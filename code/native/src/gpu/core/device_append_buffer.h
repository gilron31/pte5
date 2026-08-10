// An append-only device output buffer with an atomic bump allocator.
//
// Includable from both .cc and .cu: push() is guarded by __CUDACC__, so a host
// translation unit sees the view as a plain aggregate and is never dragged into
// nvcc. That is what keeps meet.cc ordinary C++ while push() stays a member,
// where it belongs, rather than a free function taking the view as its first
// argument.
//
// The guard does mean the class has a different definition in host and device
// translation units, which is formally an ODR violation. It is benign here and
// bounded: the members are identical, so layout and the by-value kernel
// argument are unaffected; push() is inline and __device__, so it is only ever
// emitted into device code; and there is no virtual function or static member
// for the two definitions to disagree about. This is the standard CUDA idiom --
// Thrust and CUB are built on it. Keep it that way: anything added inside the
// guard must not touch layout.
//
// This replaces passing (pointer, capacity, counter) around as three loose
// kernel arguments. The overflow contract is the part worth having in one place
// and testing directly:
//
//   *every* append bumps the counter, including the ones that do not fit.
//
// So `total()` is the number of items the kernel *wanted* to write, which is
// exactly what a caller needs to size a retry, while only the first `capacity`
// of them were stored. A caller that ignores the distinction and reads
// `capacity` items is still safe; a caller that wants completeness must compare
// total() against capacity() and retry.
#pragma once

#include "gpu/core/device_buffer.h"
#include "gpu/core/types.h"

#include <cstddef>

namespace pte::gpu {

// The device-side handle. Trivially copyable, so it is passed to kernels by
// value.
template <typename T> struct DeviceAppendView {
  T *data = nullptr;
  index_t capacity = 0;
  count_t *counter = nullptr;

#ifdef __CUDACC__
  // Appends `value`, if it fits. The counter is bumped either way; see the
  // overflow contract above.
  __device__ void push(const T &value) const {
    const count_t slot = atomicAdd(counter, count_t{1});
    if (slot < static_cast<count_t>(capacity))
      data[slot] = value;
  }
#endif
};

// The host-side owner.
template <typename T> class DeviceAppendBuffer {
public:
  explicit DeviceAppendBuffer(index_t capacity)
      : data_(static_cast<std::size_t>(capacity < 0 ? 0 : capacity)),
        counter_(1), capacity_(capacity < 0 ? 0 : capacity) {
    counter_.Zero();
    // The tail beyond what the kernel writes is never read by a correct
    // caller, but zeroing keeps a truncated copy-back deterministic and keeps
    // compute-sanitizer's initcheck quiet.
    data_.Zero();
  }

  DeviceAppendView<T> view() {
    return DeviceAppendView<T>{data_.get(), capacity_, counter_.get()};
  }

  index_t capacity() const { return capacity_; }

  // Empties the buffer for the next launch, so a caller that launches in a loop
  // (S-E runs one launch per value window) reuses one allocation instead of
  // constructing a buffer per launch -- the same argument PinnedPool makes for
  // staging memory.
  //
  // Only the counter is zeroed. The data beyond the next launch's count is
  // never read by a correct caller, and rezeroing megabytes per window to keep
  // compute-sanitizer's initcheck quiet would cost more than the launch.
  void Reset() { counter_.Zero(); }

  // Number of appends attempted, which may exceed capacity(). Synchronizes, so
  // callers that also want the data should keep the result rather than asking
  // twice.
  count_t total() const {
    count_t found = 0;
    counter_.CopyToHost(&found);
    return found;
  }

  // Number of items actually stored, given a previously read total().
  index_t stored(count_t total) const {
    return static_cast<count_t>(capacity_) < total
               ? capacity_
               : static_cast<index_t>(total);
  }

  // Copies back the first `count` items; `count` must be <= stored(total()).
  void CopyToHost(T *dst, index_t count) const {
    if (count > 0)
      data_.CopyToHost(dst, static_cast<std::size_t>(count));
  }

private:
  DeviceBuffer<T> data_;
  DeviceBuffer<count_t> counter_;
  index_t capacity_ = 0;
};

} // namespace pte::gpu
