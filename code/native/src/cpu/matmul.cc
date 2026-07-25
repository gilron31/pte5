#include "matmul.h"

#include <stdexcept>
#include <string>

namespace pte::cpu {

void matmul(const double *A, const double *B, double *C, int M, int N, int K) {
  // i-k-j order: B and C are both walked contiguously in the inner loop, which
  // is a large constant-factor win over the textbook i-j-k without changing the
  // summation order per output element.
  for (int i = 0; i < M; ++i) {
    double *c_row = C + static_cast<std::size_t>(i) * N;
    for (int j = 0; j < N; ++j)
      c_row[j] = 0.0;
    for (int k = 0; k < K; ++k) {
      const double a = A[static_cast<std::size_t>(i) * K + k];
      const double *b_row = B + static_cast<std::size_t>(k) * N;
      for (int j = 0; j < N; ++j)
        c_row[j] += a * b_row[j];
    }
  }
}

std::vector<double> matmul(const std::vector<double> &A,
                           const std::vector<double> &B, int M, int N, int K) {
  if (M < 0 || N < 0 || K < 0) {
    throw std::invalid_argument("matmul: dimensions must be non-negative");
  }
  const auto need_a = static_cast<std::size_t>(M) * K;
  const auto need_b = static_cast<std::size_t>(K) * N;
  if (A.size() != need_a || B.size() != need_b) {
    throw std::invalid_argument(
        "matmul: expected A of " + std::to_string(need_a) + " and B of " +
        std::to_string(need_b) + " elements, got " + std::to_string(A.size()) +
        " and " + std::to_string(B.size()));
  }
  std::vector<double> C(static_cast<std::size_t>(M) * N);
  matmul(A.data(), B.data(), C.data(), M, N, K);
  return C;
}

} // namespace pte::cpu
