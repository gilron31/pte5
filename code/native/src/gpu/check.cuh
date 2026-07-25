// CUDA error checking and a small RAII device buffer.
//
// Only included from .cu translation units.
#ifndef PTE_GPU_CHECK_CUH_
#define PTE_GPU_CHECK_CUH_

#include <cuda_runtime.h>

#include <cstddef>
#include <stdexcept>
#include <string>

namespace pte::gpu {

inline void ThrowIfCudaError(cudaError_t err, const char *expr,
                             const char *file, int line) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             ": " + expr + " failed: " + cudaGetErrorName(err) +
                             " (" + cudaGetErrorString(err) + ")");
  }
}

// Wrap every runtime call. Kernel launches are asynchronous, so after a launch
// use PTE_CUDA_CHECK(cudaGetLastError()) for launch-config errors and rely on
// the following synchronizing copy (also checked) for execution errors.
#define PTE_CUDA_CHECK(expr)                                                   \
  ::pte::gpu::ThrowIfCudaError((expr), #expr, __FILE__, __LINE__)

// Owning device allocation; frees on scope exit so a thrown error mid-sequence
// does not leak. cudaFree in the destructor is deliberately unchecked, since
// throwing from a destructor during unwinding would terminate.
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

  T *get() const { return ptr_; }
  std::size_t count() const { return count_; }
  std::size_t bytes() const { return count_ * sizeof(T); }

  void CopyFromHost(const T *src) {
    if (count_ > 0)
      PTE_CUDA_CHECK(cudaMemcpy(ptr_, src, bytes(), cudaMemcpyHostToDevice));
  }
  void CopyToHost(T *dst) const {
    if (count_ > 0)
      PTE_CUDA_CHECK(cudaMemcpy(dst, ptr_, bytes(), cudaMemcpyDeviceToHost));
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

#endif // PTE_GPU_CHECK_CUH_
