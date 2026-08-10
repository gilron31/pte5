// B6 -- the factor base in steady state: Map vs Load.
//
// drafts/implementation_design.md 9.4 measures construction cost and resident
// set, and explicitly makes **no throughput claim** for Map over Load because
// wall clock on this machine swings 2x. That leaves the question that actually
// decides the design bound unanswered: once the base is in use, does mapping
// cost anything per access?
//
// It might. A private file mapping and an anonymous allocation differ in page
// provenance, in whether the kernel can use transparent huge pages, and in what
// happens under memory pressure. Cycles and cache misses per norm are immune to
// the frequency scaling that blocked the earlier claim, so this is answerable
// where the earlier question was not.
//
// The access pattern is the real one -- Walker over segments -- but **strided**:
// every k-th segment, so the walk touches the whole factor base (which is the
// point, since the base is 34.6 MB at B = 1e10 against a 16 MB L3) while doing
// a fraction of the work. A contiguous prefix would measure only the first few
// megabytes and conclude, wrongly, that everything fits in cache.
//
//   ./bench_factorbase <bound> [stride]
#include "bench/harness.h"
#include "cpu/core/precision.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/generate/generate.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace pte::cpu;
using pte::bench::PerfCounters;
using pte::bench::Stopwatch;
using pte::bench::Suite;

namespace {

const char *kArtifact = "../../data/factorbase_1e12.bin";

std::size_t ResidentKiB() {
  std::FILE *f = std::fopen("/proc/self/statm", "r");
  if (!f)
    return 0;
  long total = 0, resident = 0;
  const int got = std::fscanf(f, "%ld %ld", &total, &resident);
  std::fclose(f);
  return got == 2 ? static_cast<std::size_t>(resident) * 4 : 0; // 4 KiB pages
}

struct Walked {
  std::size_t norms = 0;
  std::size_t points = 0;
};

Walked WalkStrided(const Generator &g, std::size_t stride) {
  Walked w;
  Walker walker(g);
  for (std::size_t s = 0; s < g.segment_count(); s += stride)
    walker.walk(s, [&](const NormView &n) {
      ++w.norms;
      w.points += n.points.size();
    });
  return w;
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 10'000'000'000LL;
  const std::size_t stride = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 17;

  const char *path = std::getenv("PTE_FACTORBASE");
  const std::string artifact = path ? path : kArtifact;
  const auto needed = static_cast<factor_t>(bound / kSmallestCofactor);

  std::printf("=== B6 factor base in steady state, B = %lld ===\n",
              static_cast<long long>(bound));
  std::printf("  artifact %s, primes to %llu, every %zu-th segment\n",
              artifact.c_str(), static_cast<unsigned long long>(needed),
              stride);

  const std::size_t rss0 = ResidentKiB();

  FactorBase mapped, loaded;
  double map_s = 0, load_s = 0;
  std::size_t rss_map = 0, rss_load = 0;
  try {
    Stopwatch sw;
    sw.Start();
    mapped = FactorBase::Map(artifact, needed);
    map_s = sw.WallSeconds();
    rss_map = ResidentKiB();

    sw.Start();
    loaded = FactorBase::Load(artifact, needed);
    load_s = sw.WallSeconds();
    rss_load = ResidentKiB();
  } catch (const std::exception &e) {
    std::printf("  cannot open the artifact (%s)\n", e.what());
    std::printf("  set PTE_FACTORBASE, or build one:\n"
                "    enum_campaign --bound=%lld  (it sieves its own)\n",
                static_cast<long long>(bound));
    return 0;
  }

  const double mb = static_cast<double>(mapped.size()) * 8.0 / (1024 * 1024);
  std::printf("  %zu primes, %.1f MB of records (L3 is %d MB)\n", mapped.size(),
              mb, pte::bench::Machine::kL3Bytes / (1024 * 1024));
  std::printf("  construction:  Map %.3f s (+%zu KiB RSS)   "
              "Load %.3f s (+%zu KiB RSS)\n",
              map_s, rss_map - rss0, load_s, rss_load - rss_map);

  const Generator g_map(bound, std::move(mapped));
  const Generator g_load(bound, std::move(loaded));

  // Sanity: the two bases must produce identical work, or the comparison is
  // between two different walks.
  const Walked a = WalkStrided(g_map, stride);
  const Walked b = WalkStrided(g_load, stride);
  if (a.norms != b.norms || a.points != b.points) {
    std::fprintf(stderr, "MISMATCH: mapped walk %zu norms, loaded %zu\n",
                 a.norms, b.norms);
    return 1;
  }
  std::printf("  strided walk: %zu norms, %zu points (identical both ways)\n",
              a.norms, a.points);

  Suite s("B6 walk over a mapped vs loaded factor base",
          static_cast<double>(a.norms), "norm");
  s.Note("  cycles and misses per norm; these are frequency-invariant, which "
         "is why 9.4's\n  wall-clock question is answerable here and was not "
         "there");
  s.Add("Map (mmap, file-backed)", [&] { WalkStrided(g_map, stride); });
  s.Add("Load (read into the heap)", [&] { WalkStrided(g_load, stride); });
  const auto &r = s.Run(5);

  if (r[0].counters.available) {
    const double dc =
        (static_cast<double>(r[0].counters.cycles) -
         static_cast<double>(r[1].counters.cycles)) /
        static_cast<double>(r[1].counters.cycles);
    std::printf("\n  Map costs %+.1f%% cycles against Load in steady state.\n",
                100.0 * dc);
    std::printf("  LL cache misses/norm: Map %.4f, Load %.4f\n",
                static_cast<double>(r[0].counters.cache_misses) /
                    static_cast<double>(a.norms),
                static_cast<double>(r[1].counters.cache_misses) /
                    static_cast<double>(a.norms));
  }
  return 0;
}
