// Owning device allocation.
//
// Only included from .cu translation units.
#pragma once

#include "gpu/core/check.cuh"

#include <cstddef>
#include <utility>

// GILRO: don't see a reason for it to be a .cuh and not a .h file.
namespace pte::gpu {

// Frees on scope exit, so a thrown error mid-sequence does not leak. cudaFree
// in the destructor is deliberately unchecked: throwing from a destructor
// during unwinding would terminate.
template <typename T> class DeviceBuffer {
public:
  explicit DeviceBuffer(std::size_t count) : count_(count) {
    if (count_ > 0) {
      PTE_CUDA_CHECK(cudaMalloc(&ptr_, count_ * sizeof(T)));
    }
  }
  ~DeviceBuffer() {
    if (ptr_ != nullptr)
      cudaFree(ptr_);
  }

  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;

  DeviceBuffer(DeviceBuffer &&other) noexcept
      : ptr_(std::exchange(other.ptr_, nullptr)),
        count_(std::exchange(other.count_, 0)) {}
  DeviceBuffer &operator=(DeviceBuffer &&other) noexcept {
    if (this != &other) {
      if (ptr_ != nullptr)
        cudaFree(ptr_);
      ptr_ = std::exchange(other.ptr_, nullptr);
      count_ = std::exchange(other.count_, 0);
    }
    return *this;
  }

  T *get() const { return ptr_; }
  std::size_t count() const { return count_; }
  std::size_t bytes() const { return count_ * sizeof(T); }

  void CopyFromHost(const T *src) { CopyFromHost(src, count_); }
  void CopyFromHost(const T *src, std::size_t count) {
    if (count > 0)
      PTE_CUDA_CHECK(
          cudaMemcpy(ptr_, src, count * sizeof(T), cudaMemcpyHostToDevice));
  }

  void CopyToHost(T *dst) const { CopyToHost(dst, count_); }
  // Partial copy: transferring only the prefix that was actually written keeps
  // the cost proportional to the result rather than to the capacity.
  void CopyToHost(T *dst, std::size_t count) const {
    if (count > 0)
      PTE_CUDA_CHECK(
          cudaMemcpy(dst, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost));
  }

  void Zero() {
    if (count_ > 0)
      PTE_CUDA_CHECK(cudaMemset(ptr_, 0, bytes()));
  }

private:
  T *ptr_ = nullptr;
  std::size_t count_ = 0;
};

} // namespace pte::gpu
