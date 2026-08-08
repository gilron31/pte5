// An append-only device output buffer with an atomic bump allocator.
//
// Only included from .cu translation units.
//
// This replaces passing (pointer, capacity, counter) around as three loose
// kernel arguments. The overflow contract is the part worth having in one
// place and testing directly:
//
//   *every* append bumps the counter, including the ones that do not fit.
//
// So `total()` is the number of items the kernel *wanted* to write, which is
// exactly what a caller needs to size a retry, while only the first `capacity`
// of them were stored. A caller that ignores the distinction and reads
// `capacity` items is still safe; a caller that wants completeness must
// compare total() against capacity() and retry.
#pragma once

#include "gpu/core/device_buffer.cuh"
#include "gpu/core/types.h"

#include <cstddef>

namespace pte::gpu {

// GILRO: as a rule, we should separate the "view-ish" behaviour from the
// "kernel" behaviour. They should live on two different files. If we don't it
// destroys the separation we are trying to do and enforoces unpleasant .cuh
// files everywhere where they coudl have been .h files. Let's discuss what
// should be the proper separation here.

// The device-side handle. Trivially copyable, so it is passed to kernels by
// value.

template <typename T> struct DeviceAppendView {
  T *data = nullptr;
  index_t capacity = 0;
  count_t *counter = nullptr;

  // Could add __host__
  __device__ void push(const T &value) const {
    const count_t slot = atomicAdd(counter, count_t{1});
    if (slot < static_cast<count_t>(capacity))
      data[slot] = value;
  }
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
