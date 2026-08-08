// Dense matrix multiply on the GPU. Plain C++: no <cuda_runtime.h>, no CUDA
// types.
#pragma once

namespace pte::gpu {

// Row-major (M x K) * (K x N) = (M x N), double precision, tiled kernel.
// C needs M*N doubles.
void matmul(const double *A, const double *B, double *C, int M, int N, int K);

} // namespace pte::gpu
