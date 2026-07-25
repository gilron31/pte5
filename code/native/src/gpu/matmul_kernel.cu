// Tiled shared-memory matrix multiply on the GPU, double precision, row-major.
//
// Every global read is bounds-guarded and out-of-range tile entries are zeroed,
// so dimensions that are not multiples of kTile are handled correctly.
#include "check.cuh"
#include "gpu.h"

#include <stdexcept>
#include <string>

namespace pte::gpu {
namespace {

constexpr int kTile = 16;

__global__ void MatmulKernel(const double *A, const double *B, double *C, int M,
                             int N, int K) {
  __shared__ double a_tile[kTile][kTile];
  __shared__ double b_tile[kTile][kTile];

  const int row = blockIdx.y * kTile + threadIdx.y;
  const int col = blockIdx.x * kTile + threadIdx.x;

  double sum = 0.0;
  const int tiles = (K + kTile - 1) / kTile;
  for (int t = 0; t < tiles; ++t) {
    const int a_col = t * kTile + threadIdx.x;
    const int b_row = t * kTile + threadIdx.y;

    a_tile[threadIdx.y][threadIdx.x] =
        (row < M && a_col < K) ? A[static_cast<std::size_t>(row) * K + a_col]
                               : 0.0;
    b_tile[threadIdx.y][threadIdx.x] =
        (b_row < K && col < N) ? B[static_cast<std::size_t>(b_row) * N + col]
                               : 0.0;
    __syncthreads();

    for (int k = 0; k < kTile; ++k)
      sum += a_tile[threadIdx.y][k] * b_tile[k][threadIdx.x];

    // Every thread must clear the tile before the next iteration overwrites it,
    // including threads whose own output is out of range -- hence no early
    // return anywhere in this loop.
    __syncthreads();
  }

  if (row < M && col < N)
    C[static_cast<std::size_t>(row) * N + col] = sum;
}

} // namespace

void matmul(const double *A, const double *B, double *C, int M, int N, int K) {
  if (M < 0 || N < 0 || K < 0) {
    throw std::invalid_argument("matmul: dimensions must be non-negative");
  }
  if (M == 0 || N == 0)
    return;

  const std::size_t a_count = static_cast<std::size_t>(M) * K;
  const std::size_t b_count = static_cast<std::size_t>(K) * N;
  const std::size_t c_count = static_cast<std::size_t>(M) * N;

  DeviceBuffer<double> d_a(a_count);
  DeviceBuffer<double> d_b(b_count);
  DeviceBuffer<double> d_c(c_count);
  d_a.CopyFromHost(A);
  d_b.CopyFromHost(B);

  // K == 0 means an empty sum: C is all zeros and the kernel would do nothing.
  if (K == 0) {
    d_c.Zero();
  } else {
    const dim3 threads(kTile, kTile);
    const dim3 blocks((N + kTile - 1) / kTile, (M + kTile - 1) / kTile);
    MatmulKernel<<<blocks, threads>>>(d_a.get(), d_b.get(), d_c.get(), M, N, K);
    PTE_CUDA_CHECK(cudaGetLastError());
  }
  d_c.CopyToHost(C); // synchronizes, and surfaces execution errors
}

} // namespace pte::gpu
