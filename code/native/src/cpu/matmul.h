// Dense double-precision matrix multiply on the CPU.
//
// All matrices are row-major and dimensions are (M x K) * (K x N) = (M x N),
// matching the layout expected by pte::gpu::matmul.
#ifndef PTE_CPU_MATMUL_H_
#define PTE_CPU_MATMUL_H_

#include <cstddef>
#include <vector>

namespace pte::cpu {

// C must have room for M*N doubles.
void matmul(const double *A, const double *B, double *C, int M, int N, int K);

// Convenience overload; throws std::invalid_argument on a size mismatch.
std::vector<double> matmul(const std::vector<double> &A,
                           const std::vector<double> &B, int M, int N, int K);

} // namespace pte::cpu

#endif // PTE_CPU_MATMUL_H_
