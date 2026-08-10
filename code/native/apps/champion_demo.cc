// End-to-end demo of S-E, the champion meet (drafts/gpu_implementation_plan.md
// section 3): build one enormous point set, look for a collision inside it on
// the device, and check the answer against the CPU reference.
//
// Compiled by g++, not nvcc: it reaches the kernel only through the plain-C++
// gpu/champion/champion.h.
//
// Registered with ctest as `champion_demo_cross_check`, so the CPU/GPU
// agreement is a fixture rather than something somebody has to remember to run.
// Exits non-zero on any disagreement or on a hit that fails the exact re-check.
//
//   ./champion_demo [num_primes]        default 14 (D = 2^13)
//
// The default is sized to run in a couple of seconds. The real thing is 20
// primes, D = 2^19, L1 ~ 2.5e36 and 1.374e11 pairs -- print the projection at
// the end and let the operator decide, rather than defaulting to a two-minute
// ctest case.
#include "cpu/bigint/bigint.h"
#include "cpu/champion/champion.h"
#include "cpu/champion/champion_points.h"
#include "gpu/champion/champion.h"
#include "gpu/core/device.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <span>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using pte::core::PairCollision;

double SecondsSince(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// C(D,2) as a double: the whole search's size, and the number every rate below
// is divided into.
double PairCount(std::size_t d) {
  return static_cast<double>(d) * static_cast<double>(d - 1) / 2.0;
}

std::string Decimal(const mpz_class &v) { return v.get_str(); }

void PrintHits(const std::vector<PairCollision> &hits,
               const std::vector<pte::core::U128> &nu) {
  const std::size_t show = std::min<std::size_t>(hits.size(), 4);
  for (std::size_t i = 0; i < show; ++i) {
    const PairCollision &h = hits[i];
    std::printf("    {%u,%u} | {%u,%u}   nu = %s %s %s %s\n", h.a, h.b, h.c,
                h.d, Decimal(pte::cpu::from_uint<4>(nu[h.a])).c_str(),
                Decimal(pte::cpu::from_uint<4>(nu[h.b])).c_str(),
                Decimal(pte::cpu::from_uint<4>(nu[h.c])).c_str(),
                Decimal(pte::cpu::from_uint<4>(nu[h.d])).c_str());
  }
  if (hits.size() > show)
    std::printf("    ... and %zu more\n", hits.size() - show);
}

} // namespace

int main(int argc, char **argv) {
  int num_primes = 14;
  if (argc > 1)
    num_primes = std::atoi(argv[1]);
  if (num_primes < 4 || num_primes > 21) {
    std::fprintf(stderr, "num_primes must be in [4, 21]\n");
    return 2;
  }

  const int devices = pte::gpu::device_count();
  std::printf("CUDA devices: %d\n", devices);
  if (devices > 0)
    std::printf("Using device 0: %s\n", pte::gpu::device_name(0).c_str());

  // --------------------------------------------------------------------
  // The champion norm
  // --------------------------------------------------------------------
  const auto primes = pte::cpu::smallest_split_primes(num_primes);
  const Clock::time_point t_build = Clock::now();
  const pte::cpu::ChampionNorm champ = pte::cpu::champion_norm(primes);
  const double build_s = SecondsSince(t_build);

  std::printf("\n== S-E champion meet, %d split primes ==\n", num_primes);
  std::printf("  L1                     %s\n", Decimal(champ.l1).c_str());
  std::printf("  L1 bits                %zu\n",
              mpz_sizeinbase(champ.l1.get_mpz_t(), 2));
  std::printf("  points D               %zu\n", champ.nu.size());
  std::printf("  pairs C(D,2)           %.4g\n", PairCount(champ.nu.size()));
  int nu_bits = 0;
  for (const pte::core::U128 &v : champ.nu)
    nu_bits = std::max(nu_bits, pte::core::bit_length(v));
  std::printf("  widest nu / key        %d / %d bits\n", nu_bits,
              2 * nu_bits + 1);
  std::printf("  point-set construction %.3f s   (GMP; the int32 path would "
              "overflow above 20-bit coordinates)\n",
              build_s);

  // --------------------------------------------------------------------
  // CPU reference, on a prefix
  // --------------------------------------------------------------------
  //
  // The reference is quadratic at ~1e7 pair/s, so running it on the sized
  // champion would take three hours against the device's two minutes. It is run
  // on the first `kCrossCheckPoints` points instead -- a genuine sub-point-set,
  // where both implementations must agree exactly -- while the device gets the
  // whole thing. Cross-checking a prefix is the honest trade: it checks the
  // *implementations* against each other, which is what a cross-check can do,
  // and it does not pretend the reference could have run the real search.
  constexpr std::size_t kCrossCheckPoints = 8192;
  const std::size_t cross_d = std::min(champ.nu.size(), kCrossCheckPoints);
  const std::span<const pte::core::U128> cross =
      std::span<const pte::core::U128>(champ.nu).first(cross_d);

  pte::cpu::ChampionStats cpu_stats;
  const Clock::time_point t_cpu = Clock::now();
  const std::vector<PairCollision> cpu_hits =
      pte::cpu::champion_meet(cross, {}, &cpu_stats);
  const double cpu_s = SecondsSince(t_cpu);
  std::printf("\n  CPU reference (sort each window, exact 256-bit keys), "
              "D = %zu\n",
              cross_d);
  std::printf("    windows %llu, pairs %llu, hits %llu\n",
              static_cast<unsigned long long>(cpu_stats.windows),
              static_cast<unsigned long long>(cpu_stats.pairs),
              static_cast<unsigned long long>(cpu_stats.hits));
  std::printf("    %.3f s   %.3g pair/s\n", cpu_s, PairCount(cross_d) / cpu_s);

  if (cpu_stats.pairs != static_cast<std::uint64_t>(PairCount(cross_d))) {
    std::fprintf(stderr,
                 "FAIL: the CPU reference enumerated %llu pairs, "
                 "expected C(D,2) = %.0f\n",
                 static_cast<unsigned long long>(cpu_stats.pairs),
                 PairCount(cross_d));
    return 1;
  }

  if (devices == 0) {
    std::printf("\nNo CUDA device: skipping the device half.\n");
    return 0;
  }

  // --------------------------------------------------------------------
  // Device
  // --------------------------------------------------------------------
  // Size the table to the problem rather than taking the 512 MB default: a
  // small demo should still exercise many windows instead of one launch into a
  // mostly empty table, and a full-size champion should not be cut into 65,000
  // windows, because each window costs an O(D) planning pass on the host.
  // Aim for ~32 windows, capped at the default 2^26 slots (512 MB).
  pte::gpu::ChampionParams params;
  {
    const double per_window = std::max(1.0, PairCount(champ.nu.size()) / 32.0);
    int table_log = 12;
    while (table_log < 26 &&
           static_cast<double>(1ull << table_log) * params.max_load <
               per_window)
      ++table_log;
    params.table_log = table_log;
  }
  std::printf("\n  hash table 2^%d slots (%.0f MB), window cap %.3g pairs\n",
              params.table_log,
              static_cast<double>(1ull << params.table_log) * 8 / 1048576.0,
              static_cast<double>(1ull << params.table_log) * params.max_load);
  pte::gpu::ChampionStats gpu_stats;
  // One untimed run first: the CUDA context costs ~60 ms to create and would
  // otherwise be charged to the kernel.
  (void)pte::gpu::champion_meet(
      std::span<const pte::core::U128>(champ.nu).first(8), params, nullptr);
  const Clock::time_point t_gpu = Clock::now();
  const std::vector<PairCollision> gpu_hits = pte::gpu::champion_meet(
      std::span<const pte::core::U128>(champ.nu), params, &gpu_stats);
  const double gpu_s = SecondsSince(t_gpu);

  std::printf("\n  GPU (atomicCAS hash insert per window, host confirms)\n");
  std::printf("    windows %llu, pairs %llu, candidates %llu, false positives "
              "%llu, hits %llu\n",
              static_cast<unsigned long long>(gpu_stats.windows),
              static_cast<unsigned long long>(gpu_stats.pairs),
              static_cast<unsigned long long>(gpu_stats.candidates),
              static_cast<unsigned long long>(gpu_stats.false_positives),
              static_cast<unsigned long long>(gpu_stats.hits));
  std::printf("    %.3f s   %.3g pair/s   (includes window planning on the "
              "host and %llu H2D copies)\n",
              gpu_s, PairCount(champ.nu.size()) / gpu_s,
              2ull * gpu_stats.windows);
  if (gpu_stats.pairs !=
      static_cast<std::uint64_t>(PairCount(champ.nu.size()))) {
    std::fprintf(stderr,
                 "FAIL: the device enumerated %llu pairs, expected "
                 "C(D,2) = %.0f\n",
                 static_cast<unsigned long long>(gpu_stats.pairs),
                 PairCount(champ.nu.size()));
    return 1;
  }

  // --------------------------------------------------------------------
  // The cross-check, which is the point of the program
  // --------------------------------------------------------------------
  bool ok = true;
  pte::gpu::ChampionStats gpu_cross_stats;
  const Clock::time_point t_gpu_cross = Clock::now();
  const std::vector<PairCollision> gpu_cross_hits =
      pte::gpu::champion_meet(cross, params, &gpu_cross_stats);
  const double gpu_cross_s = SecondsSince(t_gpu_cross);
  std::printf("\n  cross-check at D = %zu: CPU %.3f s, GPU %.3f s (%.0fx)\n",
              cross_d, cpu_s, gpu_cross_s, cpu_s / gpu_cross_s);

  if (gpu_cross_stats.pairs != cpu_stats.pairs) {
    std::fprintf(stderr, "FAIL: pair counts differ (CPU %llu, GPU %llu)\n",
                 static_cast<unsigned long long>(cpu_stats.pairs),
                 static_cast<unsigned long long>(gpu_cross_stats.pairs));
    ok = false;
  }
  if (cpu_hits != gpu_cross_hits) {
    std::fprintf(stderr, "FAIL: hit sets differ (CPU %zu, GPU %zu)\n",
                 cpu_hits.size(), gpu_cross_hits.size());
    ok = false;
  } else {
    std::printf("  CPU and GPU hit sets are identical -> OK  (%zu hits)\n",
                cpu_hits.size());
    PrintHits(cpu_hits, champ.nu);
  }

  // Every hit re-checked here as well, independently of either search: a
  // cross-check between two implementations proves they agree, not that they
  // are right.
  for (const PairCollision &h : gpu_hits) {
    const auto lhs = pte::core::add(pte::core::square(champ.nu[h.a]),
                                    pte::core::square(champ.nu[h.b]));
    const auto rhs = pte::core::add(pte::core::square(champ.nu[h.c]),
                                    pte::core::square(champ.nu[h.d]));
    if (!(lhs == rhs)) {
      std::fprintf(stderr, "FAIL: reported hit is not a collision\n");
      ok = false;
    }
  }

  // --------------------------------------------------------------------
  // The planted round
  // --------------------------------------------------------------------
  //
  // A champion this size has no collision, and "both sides found nothing" is a
  // fixture two broken implementations would also pass. So the same point set
  // is run again with a known collision injected: 65 = 1^2 + 8^2 = 4^2 + 7^2,
  // scaled above every real nu so the four values stay distinct from the rest.
  {
    std::vector<pte::core::U128> planted(cross.begin(), cross.end());
    const int shift = nu_bits + 2;
    const std::uint64_t base[4] = {1, 8, 4, 7};
    const std::size_t where[4] = {0, planted.size() / 3, planted.size() / 2,
                                  planted.size() - 1};
    for (int k = 0; k < 4; ++k)
      planted[where[k]] = pte::core::shl(pte::core::u128_from(base[k]), shift);

    const auto cpu_planted =
        pte::cpu::champion_meet(std::span<const pte::core::U128>(planted));
    const auto gpu_planted = pte::gpu::champion_meet(
        std::span<const pte::core::U128>(planted), params);
    const PairCollision want =
        PairCollision{static_cast<std::uint32_t>(where[0]),
                      static_cast<std::uint32_t>(where[1]),
                      static_cast<std::uint32_t>(where[2]),
                      static_cast<std::uint32_t>(where[3])}
            .canonical();
    const bool cpu_found = std::find(cpu_planted.begin(), cpu_planted.end(),
                                     want) != cpu_planted.end();
    const bool gpu_found = std::find(gpu_planted.begin(), gpu_planted.end(),
                                     want) != gpu_planted.end();
    std::printf(
        "\n  planted collision (nu ~ 2^%d, key ~ 2^%d): CPU %s, GPU %s, "
        "sets %s\n",
        shift + 3, 2 * shift + 7, cpu_found ? "found" : "MISSED",
        gpu_found ? "found" : "MISSED",
        cpu_planted == gpu_planted ? "identical" : "DIFFER");
    if (!cpu_found || !gpu_found || cpu_planted != gpu_planted)
      ok = false;
  }

  // --------------------------------------------------------------------
  // What the full-size run would cost
  // --------------------------------------------------------------------
  const double pairs_per_s = PairCount(champ.nu.size()) / gpu_s;
  std::printf("\n  at this rate, the sized champion (20 primes, D = 2^19, "
              "1.374e11 pairs) would take %.0f s\n",
              1.374e11 / pairs_per_s);
  std::printf("  bench/gpu P4 puts the hash-insert ceiling at 1.0e9 insert/s, "
              "i.e. 138 s\n");

  std::printf(ok ? "\nAll checks passed.\n" : "\nFAILURES above.\n");
  return ok ? 0 : 1;
}
