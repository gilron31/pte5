// End-to-end demo of the E_4 enumeration (search S-A of
// drafts/gpu_algorithms.md): the CPU generates admissible norms and their point
// sets, the GPU finds the key collisions, the CPU verifies every hit exactly.
//
// Compiled by g++, not nvcc: it reaches the kernel only through the plain-C++
// gpu/meet/meet.h.
//
// Registered with ctest as `enum_demo_cross_check`: it exits non-zero on any
// CPU/GPU disagreement, so the cross-check is a fixture rather than something
// somebody has to remember to run.
//
// This is deliberately the *naive hybrid* split. The timing breakdown below is
// meant to show why: the meet moves ~8 bytes per point to do a handful of
// integer operations, so the host/device transfer dominates and the GPU sits
// idle. Fixing that means generating the norms on the device -- see the
// "Deferred" section of code/native/README.md.
//
// Exits non-zero if the CPU and GPU disagree, or if any hit fails verification.
#include "cpu/confirm/confirm.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/gaussian/gaussian.h"
#include "gpu/core/device.h"
#include "gpu/meet/meet.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using pte::cpu::Gauss;
using pte::cpu::E4Sol;

// L1 bound for the demo sweep. Small enough to run in a couple of seconds and
// well under pte::cpu::kExactKeyBound, so the meet here is exact arithmetic
// rather than a filter.
constexpr std::int64_t kBound = 2000000;

double MillisSince(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start)
      .count();
}

// One admissible norm, expanded and ready to hand to either meet.
struct Work {
  std::int64_t l1 = 0;
  std::vector<pte::cpu::factor_t> primes;
  std::vector<Gauss> points;
  std::vector<pte::cpu::norm_t> nu;
};

// A verified solution, keyed so CPU and GPU results can be compared directly.
// Point indices are meaningless across implementations only if the point sets
// differ, and they do not: both sides share pte::cpu::point_set.
using Solution = std::array<std::int64_t, 5>; // l1, then the four nu values

Solution MakeSolution(const Work &w, const E4Sol &q) {
  Solution s{w.l1, w.nu[q.a], w.nu[q.b], w.nu[q.c], w.nu[q.d]};
  std::sort(s.begin() + 1, s.end());
  return s;
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

  // Context creation costs tens of milliseconds on the first launch; pay it
  // here so it is not billed to the meet.
  {
    const std::array<pte::gpu::nu_t, 4> nu = {1, 4, 7, 8};
    const std::array<pte::gpu::index_t, 2> offsets = {0, 4};
    std::array<pte::gpu::MeetHit, 1> hit{};
    const auto start = Clock::now();
    pte::gpu::meet_batch(nu, offsets, hit);
    std::printf("Device warm-up (context creation): %.2f ms\n",
                MillisSince(start));
  }

  std::printf("\n== E_4 enumeration, admissible L1 <= %lld ==\n",
              static_cast<long long>(kBound));

  // --- CPU: generate the admissible norms ----------------------------------
  auto start = Clock::now();
  const pte::cpu::Enumerator enumerator(kBound);
  std::vector<pte::cpu::DecomposedNorm> norms;
  enumerator.for_each_norm([&](const pte::cpu::DecomposedNorm &n) { norms.push_back(n); });
  const double gen_ms = MillisSince(start);

  // --- CPU: build the point sets -------------------------------------------
  start = Clock::now();
  std::vector<Work> work;
  work.reserve(norms.size());
  std::size_t total_points = 0;
  std::size_t total_pairs = 0;
  for (const pte::cpu::DecomposedNorm &n : norms) {
    Work w;
    w.l1 = n.l1;
    w.primes = enumerator.primes_of(n);
    w.points = enumerator.point_set(n);
    w.nu = pte::cpu::nu_values(w.points);
    total_points += w.points.size();
    total_pairs += w.points.size() * (w.points.size() - 1) / 2;
    work.push_back(std::move(w));
  }
  const double points_ms = MillisSince(start);

  // --- CPU: pack into the CSR layout the kernel expects --------------------
  start = Clock::now();
  std::vector<pte::gpu::nu_t> flat_nu;
  std::vector<pte::gpu::index_t> offsets;
  flat_nu.reserve(total_points);
  offsets.reserve(work.size() + 1);
  offsets.push_back(0);
  for (const Work &w : work) {
    for (const pte::cpu::norm_t v : w.nu)
      flat_nu.push_back(static_cast<pte::gpu::nu_t>(v));
    offsets.push_back(static_cast<pte::gpu::index_t>(flat_nu.size()));
  }
  const double pack_ms = MillisSince(start);

  // --- GPU: the meet -------------------------------------------------------
  start = Clock::now();
  pte::gpu::index_t capacity = 1024;
  std::vector<pte::gpu::MeetHit> raw;
  pte::gpu::index_t found = 0;
  while (true) {
    raw.assign(static_cast<std::size_t>(capacity), pte::gpu::MeetHit{});
    found = pte::gpu::meet_batch(flat_nu, offsets, raw);
    if (found <= capacity)
      break;
    capacity = found * 2; // the count reports what a retry needs
  }
  const double gpu_ms = MillisSince(start);

  // --- CPU: verify every GPU hit exactly, and keep the primitive ones -------
  start = Clock::now();
  std::vector<Solution> gpu_solutions;
  int rejected_by_verify = 0;
  int rejected_as_duplicate = 0;
  for (pte::gpu::index_t t = 0; t < found; ++t) {
    const pte::gpu::MeetHit &hit = raw[static_cast<std::size_t>(t)];
    const Work &w = work[static_cast<std::size_t>(hit.norm)];
    const E4Sol q{hit.a, hit.b, hit.c, hit.d};
    if (!pte::cpu::confirm_candidate(w.nu, q)) {
      ++rejected_by_verify;
      continue;
    }
    if (!pte::cpu::is_primitive(w.points, w.primes, q)) {
      ++rejected_as_duplicate;
      continue;
    }
    gpu_solutions.push_back(MakeSolution(w, q));
  }
  const double verify_ms = MillisSince(start);

  // --- CPU: the same search, independently, and the certificate ------------
  //
  // Every reported solution carries a machine-checked proof: the tower peeled
  // in GMP, the degree-16 identity expanded coefficient by coefficient, and
  // Theorem 4.1's divisibilities. A solution that cannot be certified is a bug
  // in this program, so it is a failure rather than a warning.
  start = Clock::now();
  std::vector<Solution> cpu_solutions;
  std::vector<pte::cpu::Certificate> certificates;
  double certify_ms = 0;
  for (const Work &w : work) {
    for (const E4Sol &q : pte::cpu::meet_oracle(w.nu)) {
      if (!pte::cpu::is_primitive(w.points, w.primes, q))
        continue;
      cpu_solutions.push_back(MakeSolution(w, q));
      const auto cert_start = Clock::now();
      certificates.push_back(pte::cpu::certify(pte::cpu::quadruple(w.points, q)));
      certify_ms += MillisSince(cert_start);
    }
  }
  const double cpu_ms = MillisSince(start) - certify_ms;

  int uncertified = 0;
  for (const pte::cpu::Certificate &c : certificates) {
    if (!c.valid()) {
      std::printf("  UNCERTIFIED: %s\n", c.failure().c_str());
      ++uncertified;
    }
  }

  std::sort(gpu_solutions.begin(), gpu_solutions.end());
  std::sort(cpu_solutions.begin(), cpu_solutions.end());

  // --- Report ---------------------------------------------------------------
  std::printf("  admissible norms   %10zu\n", work.size());
  std::printf("  points             %10zu\n", total_points);
  std::printf("  pairs              %10zu\n", total_pairs);
  std::printf("\n");
  std::printf("  CPU norm generation %8.2f ms\n", gen_ms);
  std::printf("  CPU point sets      %8.2f ms\n", points_ms);
  std::printf("  CPU pack to CSR     %8.2f ms\n", pack_ms);
  std::printf("  GPU meet            %8.2f ms  (includes H2D+D2H copies)\n",
              gpu_ms);
  std::printf("  CPU verify hits     %8.2f ms\n", verify_ms);
  std::printf("  CPU reference meet  %8.2f ms\n", cpu_ms);
  std::printf("  CPU certificates    %8.2f ms  (GMP; %zu solutions)\n",
              certify_ms, certificates.size());
  std::printf("\n");
  std::printf("  raw GPU hits        %10lld\n", static_cast<long long>(found));
  std::printf("    rejected, key mismatch in wide arithmetic %d\n",
              rejected_by_verify);
  std::printf("    rejected, imprimitive (scaled duplicate)  %d\n",
              rejected_as_duplicate);
  std::printf("  primitive E_4 solutions: GPU %zu, CPU %zu\n",
              gpu_solutions.size(), cpu_solutions.size());

  const bool agree = gpu_solutions == cpu_solutions;
  if (!agree) {
    std::printf("\n  MISMATCH between the CPU and GPU solution sets\n");
    for (const Solution &s : gpu_solutions) {
      if (!std::binary_search(cpu_solutions.begin(), cpu_solutions.end(), s))
        std::printf("    only on GPU: L1 = %lld\n",
                    static_cast<long long>(s[0]));
    }
    for (const Solution &s : cpu_solutions) {
      if (!std::binary_search(gpu_solutions.begin(), gpu_solutions.end(), s))
        std::printf("    only on CPU: L1 = %lld\n",
                    static_cast<long long>(s[0]));
    }
  } else {
    std::printf("  CPU and GPU solution sets are identical -> OK\n");
  }

  // Sorted above for the comparison, so walk the unsorted certificates by their
  // own L1 rather than assuming the two vectors still line up.
  for (const pte::cpu::Certificate &c : certificates) {
    std::printf("    E_4 at L1 = %lld  [%s]\n",
                static_cast<long long>(c.tower.l1),
                c.valid() ? "certified" : "UNCERTIFIED");
    std::printf("      L2      = %s\n", c.tower.l2.get_str().c_str());
    std::printf("      L3      = %s\n", c.tower.l3.get_str().c_str());
    std::printf("      sqrt L4 = %s\n", c.tower.sqrt_l4.get_str().c_str());
  }

  // The point of the naive-hybrid split, made numerically explicit.
  const double bytes =
      static_cast<double>(total_points) * sizeof(std::uint64_t);
  std::printf("\n  transferred %.2f MB to do ~%zu pair comparisons\n",
              bytes / (1024 * 1024), total_pairs);
  std::printf("  arithmetic intensity ~%.2f ops/byte -- transfer bound;\n",
              static_cast<double>(total_pairs) / bytes);
  std::printf("  generating norms on the device is what fixes this.\n");

  const bool ok = agree && uncertified == 0;
  if (uncertified > 0)
    std::printf("\n  %d solution(s) failed certification\n", uncertified);
  std::printf("\n%s\n", ok ? "All checks passed." : "FAILURES detected.");
  return ok ? 0 : 1;
}
