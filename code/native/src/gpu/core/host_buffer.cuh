// Owning page-locked (pinned) host allocation.
//
// Only included from .cu translation units.
//
// Why this exists, rather than just handing cudaMemcpy a std::vector:
//
//  - A copy out of pageable memory cannot be done by DMA directly. The driver
//    stages it through an internal pinned bounce buffer, so the bytes are
//    touched twice and the transfer runs at roughly half the achievable rate.
//  - cudaMemcpyAsync only actually overlaps with compute when the host side is
//    pinned; from pageable memory it silently degrades to synchronous.
//
// Two consequences shape the interface:
//
//  - cudaHostAlloc is expensive and pins physical pages, so a HostBuffer is
//    meant to be allocated once and reused across batches, never per call.
//  - Staging through it only pays off if the producer writes *into* it rather
//    than filling a std::vector and copying. Hence data()/begin()/end() and
//    operator[], so a generation loop can emit straight into the buffer.
//
// apps/hybrid_demo.cc measures the pageable-vs-pinned difference on this
// machine.
#pragma once

#include "gpu/core/check.cuh"

#include <cstddef>
#include <utility>

namespace pte::gpu {

template <typename T> class HostBuffer {
public:
  HostBuffer() = default;
  explicit HostBuffer(std::size_t count) : count_(count) {
    if (count_ > 0) {
      PTE_CUDA_CHECK(cudaHostAlloc(reinterpret_cast<void **>(&ptr_),
                                   count_ * sizeof(T), cudaHostAllocDefault));
    }
  }
  ~HostBuffer() {
    if (ptr_ != nullptr)
      cudaFreeHost(ptr_);
  }

  HostBuffer(const HostBuffer &) = delete;
  HostBuffer &operator=(const HostBuffer &) = delete;

  HostBuffer(HostBuffer &&other) noexcept
      : ptr_(std::exchange(other.ptr_, nullptr)),
        count_(std::exchange(other.count_, 0)) {}
  HostBuffer &operator=(HostBuffer &&other) noexcept {
    if (this != &other) {
      if (ptr_ != nullptr)
        cudaFreeHost(ptr_);
      ptr_ = std::exchange(other.ptr_, nullptr);
      count_ = std::exchange(other.count_, 0);
    }
    return *this;
  }

  T *data() { return ptr_; }
  const T *data() const { return ptr_; }
  std::size_t count() const { return count_; }
  std::size_t bytes() const { return count_ * sizeof(T); }

  // Range access, so a producer can fill the pinned memory directly instead of
  // building a std::vector and copying into it.
  T *begin() { return ptr_; }
  T *end() { return ptr_ + count_; }
  const T *begin() const { return ptr_; }
  const T *end() const { return ptr_ + count_; }

  T &operator[](std::size_t i) { return ptr_[i]; }
  const T &operator[](std::size_t i) const { return ptr_[i]; }

private:
  T *ptr_ = nullptr;
  std::size_t count_ = 0;
};

} // namespace pte::gpu
