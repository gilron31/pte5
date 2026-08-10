// B5 -- thread scaling of the campaign's inner loop.
//
// Two things have never been measured. First, scaling itself: 11.0x on 16
// threads is recorded in drafts/implementation_design.md 9.3, but for the
// pre-int32 build, and nothing has re-measured it since. Second, and the reason
// this bench exists at all, the **alignas(64) false-sharing fix**. ThreadState
// objects sit adjacent in a vector at ~470 bytes each, so before the fix one
// thread's per-norm counters shared a cache line with the next thread's spine
// header. ThreadSanitizer cannot see that -- it is not a race -- and a
// single-threaded profile cannot either. This sweep is the only instrument that
// can, so it runs the same work with padded and unpadded per-thread state and
// prints both.
//
// The machine has 8 cores / 16 threads, so the interesting features are the
// knee at 8 (SMT siblings start sharing an L1 and an FPU) and whether the curve
// flattens earlier than that, which would mean memory rather than compute.
//
// Bounded at B = 1e8 by default: ~0.2 s per configuration, which keeps the full
// sweep short. Larger bounds are a deliberate choice, not a default.
//
//   ./bench_scaling <bound> <max_threads>
#include "bench/harness.h"
#include "cpu/core/precision.h"
#include "cpu/generate/generate.h"
#include "cpu/meet/meet.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

using namespace pte::cpu;
using pte::bench::PerfCounters;
using pte::bench::Stopwatch;

namespace {

// The campaign's per-thread state, in both layouts. The padded one is what
// ships; the packed one is what shipped before the fix, reconstructed here so
// the difference can be measured rather than argued about.
struct Payload {
  std::size_t norms = 0;
  std::size_t hits = 0;
  MeetScratch<trunc_key_t> scratch;
  std::vector<norm_t> nu;
  Payload() : nu(kMaxPointCount) {}
};

struct Packed : Payload {};
struct alignas(64) Padded : Payload {};

template <typename State>
double RunSweep(const Generator &g, int threads, std::size_t *norms_out,
                std::size_t *hits_out) {
  std::vector<State> state(static_cast<std::size_t>(threads));
  std::atomic<std::size_t> next{0};
  const std::size_t segments = g.segment_count();

  Stopwatch sw;
  sw.Start();
  {
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));
    for (int t = 0; t < threads; ++t) {
      pool.emplace_back([&, t] {
        State &st = state[static_cast<std::size_t>(t)];
        Walker walker(g);
        for (;;) {
          const std::size_t s = next.fetch_add(1, std::memory_order_relaxed);
          if (s >= segments)
            break;
          walker.walk(s, [&](const NormView &n) {
            ++st.norms;
            for (std::size_t k = 0; k < n.points.size(); ++k)
              st.nu[k] = n.nu[k];
            meet<trunc_key_t>(n.D(), st.nu.data(), st.scratch,
                              [&](const E4Sol &) { ++st.hits; });
          });
        }
      });
    }
    for (std::thread &th : pool)
      th.join();
  }
  const double s = sw.WallSeconds(); // wall, not CPU: this is a parallel sweep

  std::size_t n = 0, h = 0;
  for (const State &st : state) {
    n += st.norms;
    h += st.hits;
  }
  *norms_out = n;
  *hits_out = h;
  return s;
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 100'000'000;
  const int hw = static_cast<int>(std::thread::hardware_concurrency());
  const int max_threads =
      argc > 2 ? std::atoi(argv[2]) : std::max(1, hw);

  std::fprintf(stderr, "building factor base for B = %lld ...\n",
               static_cast<long long>(bound));
  const Generator g(bound);

  std::printf("=== B5 thread scaling, B = %lld, %d hardware threads ===\n",
              static_cast<long long>(bound), hw);
  std::printf("  ThreadState layouts: Packed (%zu B, pre-fix) vs "
              "alignas(64) Padded (%zu B, shipped)\n",
              sizeof(Packed), sizeof(Padded));
  std::printf("\n  %8s %12s %12s %10s %10s %9s %9s\n", "threads", "packed s",
              "padded s", "packed Mn/s", "padded Mn/s", "speedup", "padding");

  std::size_t base_norms = 0, hits = 0;
  double t1_padded = 0;
  bool consistent = true;

  for (int t = 1; t <= max_threads; t *= 2) {
    // Interleave the two layouts within a thread count, and take the best of
    // three: the layouts must see the same thermal conditions or the
    // comparison is worthless.
    double packed = 1e300, padded = 1e300;
    std::size_t n1 = 0, n2 = 0, h1 = 0, h2 = 0;
    for (int r = 0; r < 3; ++r) {
      packed = std::min(packed, RunSweep<Packed>(g, t, &n1, &h1));
      padded = std::min(padded, RunSweep<Padded>(g, t, &n2, &h2));
    }
    if (base_norms == 0) {
      base_norms = n1;
      hits = h1;
      t1_padded = padded;
    }
    // Every configuration must enumerate exactly the same norms; a scaling
    // curve over a work total that changes with thread count is meaningless.
    if (n1 != base_norms || n2 != base_norms || h1 != hits || h2 != hits)
      consistent = false;

    const double mn = static_cast<double>(base_norms) / 1e6;
    std::printf("  %8d %12.3f %12.3f %10.2f %10.2f %8.2fx %8.1f%%\n", t, packed,
                padded, mn / packed, mn / padded, t1_padded / padded,
                100.0 * (packed - padded) / packed);
  }

  if (max_threads > 1 && (max_threads & (max_threads - 1)) != 0) {
    double packed = 1e300, padded = 1e300;
    std::size_t n1 = 0, n2 = 0, h1 = 0, h2 = 0;
    for (int r = 0; r < 3; ++r) {
      packed = std::min(packed, RunSweep<Packed>(g, max_threads, &n1, &h1));
      padded = std::min(padded, RunSweep<Padded>(g, max_threads, &n2, &h2));
    }
    const double mn = static_cast<double>(base_norms) / 1e6;
    std::printf("  %8d %12.3f %12.3f %10.2f %10.2f %8.2fx %8.1f%%\n",
                max_threads, packed, padded, mn / packed, mn / padded,
                t1_padded / padded, 100.0 * (packed - padded) / packed);
    if (n1 != base_norms || n2 != base_norms)
      consistent = false;
  }

  std::printf("\n  %zu norms, %zu meet hits, identical at every thread count: "
              "%s\n",
              base_norms, hits, consistent ? "yes" : "NO -- results above are "
                                                     "not comparable");
  std::printf("  'padding' is how much the alignas(64) layout wins; a figure "
              "in the noise at 1-2\n"
              "  threads and a positive one at 8-16 is the signature of false "
              "sharing.\n");
  return consistent ? 0 : 1;
}
