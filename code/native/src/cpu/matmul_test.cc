// Unit tests for the CPU matrix multiply.
//
// The oracle is Boost.uBLAS's prod(), which shares no code with matmul.cc.
#include "matmul.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

namespace {

namespace ublas = boost::numeric::ublas;

// uBLAS matrix<double> is row-major by default, same as our convention.
ublas::matrix<double> ToUblas(const std::vector<double> &v, int rows,
                              int cols) {
  ublas::matrix<double> m(rows, cols);
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      m(i, j) = v[static_cast<std::size_t>(i) * cols + j];
  return m;
}

std::vector<double> RandomMatrix(int rows, int cols, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::vector<double> v(static_cast<std::size_t>(rows) * cols);
  for (auto &x : v)
    x = dist(rng);
  return v;
}

// Both implementations accumulate in a different order, so exact equality is
// not required -- only that the gap stays at rounding level.
void ExpectMatchesUblas(int M, int N, int K, unsigned seed) {
  const auto A = RandomMatrix(M, K, seed);
  const auto B = RandomMatrix(K, N, seed + 1);
  const auto C = pte::cpu::matmul(A, B, M, N, K);

  const ublas::matrix<double> expected =
      ublas::prod(ToUblas(A, M, K), ToUblas(B, K, N));

  ASSERT_EQ(C.size(), static_cast<std::size_t>(M) * N);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      const double got = C[static_cast<std::size_t>(i) * N + j];
      EXPECT_NEAR(got, expected(i, j), 1e-12 * (1.0 + std::abs(expected(i, j))))
          << "M=" << M << " N=" << N << " K=" << K << " at (" << i << "," << j
          << ")";
    }
  }
}

TEST(CpuMatmul, MatchesUblasSquare) { ExpectMatchesUblas(32, 32, 32, 1); }

TEST(CpuMatmul, MatchesUblasNonSquare) {
  ExpectMatchesUblas(7, 13, 5, 2);
  ExpectMatchesUblas(13, 7, 29, 3);
  ExpectMatchesUblas(1, 1, 1, 4);
  ExpectMatchesUblas(1, 64, 64, 5);  // row vector times matrix
  ExpectMatchesUblas(64, 1, 64, 6);  // matrix times column vector
  ExpectMatchesUblas(37, 53, 41, 7); // all three prime, no nice blocking
}

TEST(CpuMatmul, IdentityIsNeutral) {
  const int n = 16;
  const auto A = RandomMatrix(n, n, 42);
  std::vector<double> I(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    I[static_cast<std::size_t>(i) * n + i] = 1.0;

  const auto AI = pte::cpu::matmul(A, I, n, n, n);
  const auto IA = pte::cpu::matmul(I, A, n, n, n);
  for (std::size_t i = 0; i < A.size(); ++i) {
    EXPECT_DOUBLE_EQ(AI[i], A[i]) << "i=" << i;
    EXPECT_DOUBLE_EQ(IA[i], A[i]) << "i=" << i;
  }
}

TEST(CpuMatmul, KnownSmallProduct) {
  // [[1,2,3],[4,5,6]] * [[7,8],[9,10],[11,12]] = [[58,64],[139,154]]
  const std::vector<double> A = {1, 2, 3, 4, 5, 6};
  const std::vector<double> B = {7, 8, 9, 10, 11, 12};
  const auto C = pte::cpu::matmul(A, B, 2, 2, 3);
  ASSERT_EQ(C.size(), 4u);
  EXPECT_DOUBLE_EQ(C[0], 58.0);
  EXPECT_DOUBLE_EQ(C[1], 64.0);
  EXPECT_DOUBLE_EQ(C[2], 139.0);
  EXPECT_DOUBLE_EQ(C[3], 154.0);
}

TEST(CpuMatmul, DoesNotAssumeCIsZeroed) {
  // matmul must overwrite C, not accumulate into whatever was there.
  const std::vector<double> A = {1, 2, 3, 4};
  const std::vector<double> B = {1, 0, 0, 1};
  std::vector<double> C(4, 999.0);
  pte::cpu::matmul(A.data(), B.data(), C.data(), 2, 2, 2);
  EXPECT_DOUBLE_EQ(C[0], 1.0);
  EXPECT_DOUBLE_EQ(C[1], 2.0);
  EXPECT_DOUBLE_EQ(C[2], 3.0);
  EXPECT_DOUBLE_EQ(C[3], 4.0);
}

TEST(CpuMatmul, RejectsSizeMismatch) {
  const std::vector<double> A = {1, 2, 3, 4};
  const std::vector<double> B = {1, 2};
  EXPECT_THROW(pte::cpu::matmul(A, B, 2, 2, 2), std::invalid_argument);
}

} // namespace
