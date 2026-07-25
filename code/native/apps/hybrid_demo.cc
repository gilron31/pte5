// End-to-end hybrid demo: the CPU prepares work with GMP, the GPU does the
// arithmetic, the CPU verifies the results against GMP again.
//
// This is the program that ties the two halves together: it is compiled by g++
// (not nvcc) and reaches the kernels only through the plain-C++ src/gpu/gpu.h.
//
// Exits non-zero if any GPU result disagrees with the CPU reference.
#include "bigint.h"
#include "gpu.h"
#include "matmul.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double MillisSince(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// --- Use case 1: batched big-integer multiply -------------------------------
//
// GMP generates the operands and checks the products; the GPU multiplies them
// with the naive schoolbook kernel. The limb array is the only thing crossing
// between the two.
bool RunBigIntDemo() {
  constexpr int kBits = 256;
  constexpr int kLimbs = kBits / pte::cpu::kLimbBits; // 8
  constexpr int kBatch = 20000;

  std::printf("\n== Big integers: %d pairs of %d-bit operands ==\n", kBatch,
              kBits);

  // CPU: generate the operands with GMP.
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(20260725);
  std::vector<mpz_class> a(kBatch), b(kBatch);
  for (int t = 0; t < kBatch; ++t) {
    a[t] = rng.get_z_bits(kBits);
    b[t] = rng.get_z_bits(kBits);
  }

  // CPU reference: multiply with GMP and time it.
  auto start = Clock::now();
  std::vector<mpz_class> expected(kBatch);
  for (int t = 0; t < kBatch; ++t)
    expected[t] = pte::cpu::mul(a[t], b[t]);
  const double cpu_ms = MillisSince(start);

  // CPU: flatten into the limb layout the kernel expects.
  start = Clock::now();
  std::vector<std::uint32_t> flat_a(static_cast<std::size_t>(kBatch) * kLimbs);
  std::vector<std::uint32_t> flat_b(static_cast<std::size_t>(kBatch) * kLimbs);
  for (int t = 0; t < kBatch; ++t) {
    pte::cpu::to_limbs(a[t], flat_a.data() + static_cast<std::size_t>(t) * kLimbs,
                       kLimbs);
    pte::cpu::to_limbs(b[t], flat_b.data() + static_cast<std::size_t>(t) * kLimbs,
                       kLimbs);
  }
  const double pack_ms = MillisSince(start);

  // GPU: multiply the whole batch. Timing includes the host/device transfers,
  // which at this size dominate the arithmetic.
  start = Clock::now();
  std::vector<std::uint32_t> flat_out(static_cast<std::size_t>(kBatch) * 2 *
                                      kLimbs);
  pte::gpu::mul_batch(flat_a.data(), flat_b.data(), flat_out.data(), kLimbs,
                      kBatch);
  const double gpu_ms = MillisSince(start);

  // CPU: unpack and verify every single product against GMP.
  start = Clock::now();
  int mismatches = 0;
  for (int t = 0; t < kBatch; ++t) {
    const mpz_class got = pte::cpu::from_limbs(
        flat_out.data() + static_cast<std::size_t>(t) * 2 * kLimbs, 2 * kLimbs);
    if (got != expected[t]) {
      if (mismatches < 3) {
        std::printf("  MISMATCH at %d\n    gpu = %s\n    gmp = %s\n", t,
                    got.get_str().c_str(), expected[t].get_str().c_str());
      }
      ++mismatches;
    }
  }
  const double verify_ms = MillisSince(start);

  std::printf("  CPU (GMP mul)      %8.2f ms\n", cpu_ms);
  std::printf("  CPU pack to limbs  %8.2f ms\n", pack_ms);
  std::printf("  GPU (schoolbook)   %8.2f ms  (includes H2D+D2H copies)\n",
              gpu_ms);
  std::printf("  CPU verify vs GMP  %8.2f ms\n", verify_ms);
  std::printf("  mismatches: %d / %d -> %s\n", mismatches, kBatch,
              mismatches == 0 ? "OK" : "FAILED");
  return mismatches == 0;
}

// --- Use case 2: matrix multiply --------------------------------------------
bool RunMatmulDemo() {
  constexpr int kM = 512, kN = 512, kK = 512;
  std::printf("\n== Matrix multiply: (%dx%d) * (%dx%d), double ==\n", kM, kK, kK,
              kN);

  std::mt19937 rng(7);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::vector<double> A(static_cast<std::size_t>(kM) * kK);
  std::vector<double> B(static_cast<std::size_t>(kK) * kN);
  for (auto &x : A)
    x = dist(rng);
  for (auto &x : B)
    x = dist(rng);

  auto start = Clock::now();
  const auto cpu_c = pte::cpu::matmul(A, B, kM, kN, kK);
  const double cpu_ms = MillisSince(start);

  start = Clock::now();
  std::vector<double> gpu_c(static_cast<std::size_t>(kM) * kN);
  pte::gpu::matmul(A.data(), B.data(), gpu_c.data(), kM, kN, kK);
  const double gpu_ms = MillisSince(start);

  // Tile-order versus row-order accumulation, so expect rounding-level gaps.
  double max_abs_diff = 0.0;
  for (std::size_t i = 0; i < cpu_c.size(); ++i)
    max_abs_diff = std::max(max_abs_diff, std::abs(cpu_c[i] - gpu_c[i]));

  const double gflop = 2.0 * kM * kN * kK / 1e9;
  std::printf("  CPU (i-k-j loops)  %8.2f ms  (%.2f GFLOP/s)\n", cpu_ms,
              gflop / (cpu_ms / 1000.0));
  std::printf("  GPU (16x16 tiled)  %8.2f ms  (%.2f GFLOP/s, incl. copies)\n",
              gpu_ms, gflop / (gpu_ms / 1000.0));
  std::printf("  max |cpu - gpu|: %.3e -> %s\n", max_abs_diff,
              max_abs_diff < 1e-9 ? "OK" : "FAILED");
  return max_abs_diff < 1e-9;
}

} // namespace

int main() {
  const int devices = pte::gpu::device_count();
  std::printf("CUDA devices: %d\n", devices);
  if (devices == 0) {
    std::printf("No CUDA device available; nothing to demonstrate.\n");
    return 1;
  }
  std::printf("Using device 0: %s\n", pte::gpu::device_name(0).c_str());

  // Creating the CUDA context costs tens of milliseconds and happens on the
  // first kernel launch. Do it here so it is not billed to the first demo.
  {
    const std::uint32_t one = 1;
    std::uint32_t scratch[2] = {0, 0};
    auto start = Clock::now();
    pte::gpu::mul_batch(&one, &one, scratch, 1, 1);
    std::printf("Device warm-up (context creation): %.2f ms\n",
                MillisSince(start));
  }

  const bool bigint_ok = RunBigIntDemo();
  const bool matmul_ok = RunMatmulDemo();

  std::printf("\n%s\n", (bigint_ok && matmul_ok) ? "All checks passed."
                                                 : "FAILURES detected.");
  return (bigint_ok && matmul_ok) ? 0 : 1;
}
