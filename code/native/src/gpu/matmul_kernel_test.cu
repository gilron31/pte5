// Unit tests for the GPU tiled matrix multiply.
//
// Oracle is a host triple loop written here, so this binary depends on pte_gpu
// alone -- no GMP, no pte_cpu.
#include "gpu.h"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

namespace {

class GpuMatmul : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

std::vector<double> HostMatmul(const std::vector<double> &A,
                              const std::vector<double> &B, int M, int N,
                              int K) {
  std::vector<double> C(static_cast<std::size_t>(M) * N, 0.0);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j) {
      double sum = 0.0;
      for (int k = 0; k < K; ++k)
        sum += A[static_cast<std::size_t>(i) * K + k] *
               B[static_cast<std::size_t>(k) * N + j];
      C[static_cast<std::size_t>(i) * N + j] = sum;
    }
  return C;
}

std::vector<double> RandomMatrix(int rows, int cols, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::vector<double> v(static_cast<std::size_t>(rows) * cols);
  for (auto &x : v)
    x = dist(rng);
  return v;
}

// The kernel accumulates over tiles while the host loop runs straight through,
// so agreement is only up to rounding.
void ExpectMatchesHost(int M, int N, int K, unsigned seed) {
  const auto A = RandomMatrix(M, K, seed);
  const auto B = RandomMatrix(K, N, seed + 1);
  const auto expected = HostMatmul(A, B, M, N, K);

  std::vector<double> got(static_cast<std::size_t>(M) * N, 0.0);
  pte::gpu::matmul(A.data(), B.data(), got.data(), M, N, K);

  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j) {
      const std::size_t idx = static_cast<std::size_t>(i) * N + j;
      EXPECT_NEAR(got[idx], expected[idx],
                  1e-12 * (1.0 + std::abs(expected[idx])))
          << "M=" << M << " N=" << N << " K=" << K << " at (" << i << "," << j
          << ")";
    }
}

TEST_F(GpuMatmul, MatchesHostOnTileMultiples) {
  ExpectMatchesHost(16, 16, 16, 1);
  ExpectMatchesHost(32, 64, 48, 2);
}

// The classic tiled-kernel bug is a missing bounds guard, which only shows up
// when a dimension is not a multiple of the tile size.
TEST_F(GpuMatmul, MatchesHostOnRaggedDimensions) {
  ExpectMatchesHost(37, 53, 41, 3); // all three prime
  ExpectMatchesHost(1, 1, 1, 4);
  ExpectMatchesHost(17, 15, 33, 5);  // straddles one tile in every dimension
  ExpectMatchesHost(1, 100, 100, 6); // single row
  ExpectMatchesHost(100, 1, 100, 7); // single column
  ExpectMatchesHost(100, 100, 1, 8); // rank-one product, K < kTile
  ExpectMatchesHost(15, 15, 15, 9);  // smaller than one tile in all dimensions
}

TEST_F(GpuMatmul, MatchesHostOnBigMatrix) {
  ExpectMatchesHost(200, 173, 191, 10);
}

TEST_F(GpuMatmul, KnownSmallProduct) {
  // [[1,2,3],[4,5,6]] * [[7,8],[9,10],[11,12]] = [[58,64],[139,154]]
  const std::vector<double> A = {1, 2, 3, 4, 5, 6};
  const std::vector<double> B = {7, 8, 9, 10, 11, 12};
  std::vector<double> C(4, -1.0);
  pte::gpu::matmul(A.data(), B.data(), C.data(), 2, 2, 3);
  EXPECT_DOUBLE_EQ(C[0], 58.0);
  EXPECT_DOUBLE_EQ(C[1], 64.0);
  EXPECT_DOUBLE_EQ(C[2], 139.0);
  EXPECT_DOUBLE_EQ(C[3], 154.0);
}

TEST_F(GpuMatmul, IdentityIsNeutral) {
  const int n = 40; // not a tile multiple
  const auto A = RandomMatrix(n, n, 99);
  std::vector<double> I(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    I[static_cast<std::size_t>(i) * n + i] = 1.0;

  std::vector<double> C(static_cast<std::size_t>(n) * n, 0.0);
  pte::gpu::matmul(A.data(), I.data(), C.data(), n, n, n);
  for (std::size_t i = 0; i < A.size(); ++i)
    EXPECT_DOUBLE_EQ(C[i], A[i]) << "i=" << i;
}

TEST_F(GpuMatmul, OverwritesRatherThanAccumulates) {
  const std::vector<double> A = {1, 2, 3, 4};
  const std::vector<double> B = {1, 0, 0, 1};
  std::vector<double> C(4, 999.0);
  pte::gpu::matmul(A.data(), B.data(), C.data(), 2, 2, 2);
  EXPECT_DOUBLE_EQ(C[0], 1.0);
  EXPECT_DOUBLE_EQ(C[1], 2.0);
  EXPECT_DOUBLE_EQ(C[2], 3.0);
  EXPECT_DOUBLE_EQ(C[3], 4.0);
}

TEST_F(GpuMatmul, ZeroKGivesZeroMatrix) {
  std::vector<double> C(4, 999.0);
  pte::gpu::matmul(nullptr, nullptr, C.data(), 2, 2, 0);
  for (double x : C)
    EXPECT_DOUBLE_EQ(x, 0.0);
}

TEST(GpuMatmulArgs, RejectsNegativeDimensions) {
  // Validated before any device work, so this needs no GPU.
  std::vector<double> C(4, 0.0);
  EXPECT_THROW(pte::gpu::matmul(nullptr, nullptr, C.data(), -1, 2, 2),
               std::invalid_argument);
  EXPECT_THROW(pte::gpu::matmul(nullptr, nullptr, C.data(), 2, 2, -3),
               std::invalid_argument);
}

} // namespace
