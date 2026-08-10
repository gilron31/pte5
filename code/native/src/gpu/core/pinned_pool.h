// Pooled pinned staging: one page-locked allocation, reused for every transfer.
//
// ---------------------------------------------------------------------------
// The defect this closes
// ---------------------------------------------------------------------------
//
// `meet.cc` stages each batch through a fresh `HostBuffer`, i.e. it calls
// cudaHostAlloc once per launch. REVIEW.md flagged that as a regression when it
// landed (the demo's GPU stage went 0.92 ms -> 3.77 ms) and left it open;
// bench/gpu P7 then measured it on this machine:
//
//   cudaHostAlloc per transfer   1.45 GB/s
//   pooled pinned buffer        13.4 GB/s     <- 9.2x
//   pageable host memory        13.1 GB/s
//
// Read the third row before drawing the usual conclusion: **the win is not
// pinning, it is not allocating.** Pageable memory is within 3% of pooled
// pinned for a synchronous copy, because the driver's own bounce buffer is
// itself pooled. Allocating fresh pinned memory per transfer is nine times
// *worse* than doing nothing at all -- cudaHostAlloc pins physical pages, which
// means a kernel-side page-table walk and a TLB shootdown every call.
//
// So the reason to keep pinned memory is not this measurement. It is that
// cudaMemcpyAsync only overlaps with compute from pinned memory (from pageable
// it silently degrades to synchronous), and that is the door S-E's per-window
// pipeline wants left open. The reason to *pool* it is the 9.2x.
//
// ---------------------------------------------------------------------------
// The contract
// ---------------------------------------------------------------------------
//
// One live region at a time. `Acquire<T>(n)` hands back a span into the single
// pooled block, and the next `Acquire` invalidates it. That is deliberately the
// weakest contract that supports the intended use -- fill the staging buffer,
// issue a *synchronous* copy, repeat -- because a pool that hands out several
// independent regions is an arena, and an arena wants a lifetime story that
// nothing here needs yet. Two regions that must be live together come from one
// Acquire and a `subspan`; S-E's per-window lo/hi arrays do exactly that.
//
// The block only ever grows, so a steady-state loop allocates nothing:
// `allocations()` is exposed so a test can assert that rather than infer it.
//
// Producers should write *into* the acquired span rather than filling a vector
// and copying it in -- REVIEW.md A2(a), and the reason core/window.h's
// WindowRanges holds spans instead of vectors.
#pragma once

#include "gpu/core/device_buffer.h"
#include "gpu/core/host_buffer.h"

#include <cstddef>
#include <cstdint>
#include <span>

namespace pte::gpu {

class PinnedPool {
public:
  PinnedPool() = default;

  PinnedPool(const PinnedPool &) = delete;
  PinnedPool &operator=(const PinnedPool &) = delete;

  // A staging region of `count` T's in pinned memory. Invalidated by the next
  // Acquire.
  //
  // Alignment: cudaHostAlloc returns page-aligned memory (4096 B here), so the
  // cast is well aligned for every type this project stages -- the widest is
  // core::U256 at 8. The static_assert makes that a check rather than a
  // sentence.
  template <typename T> std::span<T> Acquire(std::size_t count) {
    static_assert(alignof(T) <= 256,
                  "PinnedPool assumes cudaHostAlloc's page alignment covers T");
    Reserve(count * sizeof(T));
    if (count == 0)
      return {};
    return std::span<T>(reinterpret_cast<T *>(block_.data()), count);
  }

  // Stage a pageable source and copy it to the device. The convenience form,
  // for callers that already have the data somewhere else; a producer that can
  // emit straight into pinned memory should use Acquire instead and save the
  // extra pass.
  template <typename T>
  void StageAndCopy(const T *src, std::size_t count, DeviceBuffer<T> &dst) {
    if (count == 0)
      return;
    const std::span<T> staging = Acquire<T>(count);
    for (std::size_t i = 0; i < count; ++i)
      staging[i] = src[i];
    dst.CopyFromHost(staging.data(), count);
  }

  // Bytes currently pinned, and the number of cudaHostAlloc calls made so far.
  // The second is the point of the class: in a loop it must stop rising.
  std::size_t bytes() const { return block_.bytes(); }
  std::size_t allocations() const { return allocations_; }

private:
  // Geometric growth, so a loop whose sizes drift upwards still converges to a
  // constant number of allocations instead of one per size.
  void Reserve(std::size_t bytes) {
    if (bytes <= block_.bytes())
      return;
    std::size_t want = block_.bytes() == 0 ? 4096 : block_.bytes();
    while (want < bytes)
      want *= 2;
    block_ = HostBuffer<unsigned char>(want);
    ++allocations_;
  }

  HostBuffer<unsigned char> block_;
  std::size_t allocations_ = 0;
};

} // namespace pte::gpu
