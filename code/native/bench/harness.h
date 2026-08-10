// The measurement discipline, in one place.
//
// This machine is a laptop whose cores range 0.4-5.1 GHz under an `amd-pstate`
// governor. Wall-clock runs of the same work have differed by 2x, and that
// variance has already forced one published claim in
// drafts/implementation_design.md to be retracted. Everything here exists to
// make a number defensible anyway:
//
//   1. INTERLEAVE.  Variants run round-robin, never one-after-the-other. A
//      thermal ramp or a governor step then hits every variant equally, so the
//      *ratio* survives even when the absolutes drift. Ratios are what
//      benchmarks are actually for.
//
//   2. REPORT THE SPREAD.  min, median and inter-quartile range. The IQR is the
//      error bar; a result whose IQR exceeds the difference it claims to show
//      is not a result.
//
//   3. WEIGH THE ENGINE.  A fixed dependent-latency chain of known cycle cost
//      runs before and after every sweep. It measures the clock the sweep
//      actually ran at, and its drift across the sweep is printed. More than a
//      few percent and the table says UNSTABLE rather than pretending.
//
//   4. COUNT, DO NOT TIME, WHERE POSSIBLE.  With perf_event_open available
//      (see perf_counters.h) cycles-per-unit replaces ns-per-unit as the
//      headline: it is invariant under frequency scaling, which removes the
//      whole problem rather than working around it.
//
// The dynamometer earns its keep even after `cpupower frequency-set -g
// performance`: it is how you find out the machine thermally throttled halfway
// through a sweep, which no amount of governor pinning prevents.
#pragma once

#include "bench/ceiling.h"
#include "bench/perf_counters.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace pte::bench {

// ---------------------------------------------------------------------------
// Time
// ---------------------------------------------------------------------------

// Wall time and this thread's CPU time. They diverge exactly when something
// else is competing for the core, which is a fact worth knowing rather than
// averaging away.
class Stopwatch {
public:
  void Start() {
    wall_ = std::chrono::steady_clock::now();
    cpu_ = ThreadCpuSeconds();
  }
  double WallSeconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         wall_)
        .count();
  }
  double CpuSeconds() const { return ThreadCpuSeconds() - cpu_; }

  static double ThreadCpuSeconds() {
    timespec ts{};
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return static_cast<double>(ts.tv_sec) +
           static_cast<double>(ts.tv_nsec) * 1e-9;
  }

private:
  std::chrono::steady_clock::time_point wall_{};
  double cpu_ = 0;
};

// ---------------------------------------------------------------------------
// The dynamometer
// ---------------------------------------------------------------------------

// A dependent chain of 64-bit multiply-add. Each iteration waits on the
// previous one, so the loop is latency-bound and its rate is proportional to
// core frequency and to nothing else -- no ILP, no cache, no branch prediction.
// That is what makes it a clock probe rather than a benchmark.
inline std::uint64_t SpinChain(std::uint64_t iters, std::uint64_t seed = 1) {
  std::uint64_t x = seed;
  for (std::uint64_t i = 0; i < iters; ++i)
    x = x * 6364136223846793005ULL + 1442695040888963407ULL;
  return x;
}

class Dynamometer {
public:
  struct Reading {
    double ghz = 0;
    bool exact = false; // true when measured from the cycle counter
  };

  // ~0.4 ms of work per sample, best of three.
  //
  // Two details that both cost a day if you get them wrong. **Thread CPU time,
  // not wall time**: the window is short enough that one scheduler tick
  // inflates wall while the cycle counter (user-space only) does not move,
  // which reads as a 30% clock drop that never happened -- measured, on this
  // machine, before the fix. And **best of three**, because the failure mode is
  // one-sided: interference can only make the apparent clock lower, never
  // higher.
  static constexpr std::uint64_t kIters = 1u << 19;
  static constexpr int kSamples = 3;

  Reading Measure() {
    Reading best;
    for (int i = 0; i < kSamples; ++i) {
      PerfCounters pc;
      Stopwatch sw;
      pc.Start();
      sw.Start();
      sink_ += SpinChain(kIters, static_cast<std::uint64_t>(i) + 1);
      const double s = sw.CpuSeconds();
      const PerfCounters::Reading r = pc.Stop();
      if (s <= 0)
        continue;
      const Reading here =
          (r.available && r.cycles > 0)
              ? Reading{static_cast<double>(r.cycles) / s / 1e9, true}
              // Fallback: imul(3) + add(1) on the critical path.
              : Reading{static_cast<double>(kIters) * 4.0 / s / 1e9, false};
      if (here.ghz > best.ghz)
        best = here;
    }
    return best;
  }

  std::uint64_t sink() const { return sink_; }

private:
  std::uint64_t sink_ = 0;
};

// ---------------------------------------------------------------------------
// The runner
// ---------------------------------------------------------------------------

struct Result {
  std::string name;
  double min_s = 0;
  double median_s = 0;
  double iqr_frac = 0; // (q3 - q1) / median
  PerfCounters::Reading counters;
};

class Suite {
public:
  // `units` is what one pass of a variant processes -- norms, partitions, keys.
  // Every reported rate is per unit, which is what makes bounds comparable.
  Suite(std::string title, double units, std::string unit_name)
      : title_(std::move(title)), units_(units),
        unit_name_(std::move(unit_name)) {}

  Suite &SetCeiling(Ceiling c) {
    ceiling_ = std::move(c);
    return *this;
  }
  // Extra free-form line under the title, e.g. corpus provenance.
  Suite &Note(std::string n) {
    notes_.push_back(std::move(n));
    return *this;
  }

  void Add(std::string name, std::function<void()> body) {
    variants_.push_back({std::move(name), std::move(body)});
  }

  // Runs every variant `repeats` times, interleaved, then prints the table.
  const std::vector<Result> &Run(int repeats) {
    const auto n = variants_.size();
    std::vector<std::vector<double>> times(n);
    std::vector<PerfCounters::Reading> best(n);
    std::vector<double> best_s(n, 1e300);

    Dynamometer dyno;

    // Warm-up round, discarded. Without it the first calibration is taken on a
    // core still at its idle P-state, and the "drift" that follows is really
    // the governor ramping -- measured at 31% on this machine, which would
    // stamp UNSTABLE on every otherwise-clean sweep. It also faults in the
    // corpus and fills the caches, so the min over the measured rounds means
    // "steady state" rather than "whichever round got lucky".
    for (const Variant &v : variants_)
      v.body();
    const Dynamometer::Reading ghz_before = dyno.Measure();
    double ghz_min = ghz_before.ghz, ghz_max = ghz_before.ghz;

    // One pass over all variants per round: a governor step or a thermal event
    // lands on every variant, not on whichever ran last.
    for (int r = 0; r < repeats; ++r) {
      for (std::size_t v = 0; v < n; ++v) {
        PerfCounters pc;
        Stopwatch sw;
        pc.Start();
        sw.Start();
        variants_[v].body();
        const double s = sw.WallSeconds();
        const PerfCounters::Reading cr = pc.Stop();
        times[v].push_back(s);
        if (s < best_s[v]) {
          best_s[v] = s;
          best[v] = cr;
        }
      }
      const Dynamometer::Reading g = dyno.Measure();
      ghz_min = std::min(ghz_min, g.ghz);
      ghz_max = std::max(ghz_max, g.ghz);
    }

    results_.clear();
    for (std::size_t v = 0; v < n; ++v) {
      std::vector<double> t = times[v];
      std::sort(t.begin(), t.end());
      const double med = t[t.size() / 2];
      const double q1 = t[t.size() / 4];
      const double q3 = t[(3 * t.size()) / 4];
      results_.push_back(Result{variants_[v].name, t.front(), med,
                                med > 0 ? (q3 - q1) / med : 0.0, best[v]});
    }

    ghz_ = ghz_max;
    ghz_exact_ = ghz_before.exact;
    drift_ = ghz_max > 0 ? (ghz_max - ghz_min) / ghz_max : 0.0;
    sink_ += dyno.sink();
    Report();
    return results_;
  }

  const std::vector<Result> &results() const { return results_; }
  double ghz() const { return ghz_; }
  std::uint64_t sink() const { return sink_; }

private:
  struct Variant {
    std::string name;
    std::function<void()> body;
  };

  void Report() const {
    std::printf("\n=== %s ===\n", title_.c_str());
    std::printf("%.0f %s per pass", units_, unit_name_.c_str());
    std::printf("   |   clock %.2f GHz%s", ghz_, ghz_exact_ ? "" : " (est.)");
    std::printf("   |   drift %.1f%% %s\n", 100.0 * drift_,
                drift_ > 0.03 ? "<-- UNSTABLE, treat absolutes with suspicion"
                              : "(stable)");
    for (const std::string &s : notes_)
      std::printf("%s\n", s.c_str());
    if (!ceiling_.empty())
      std::printf("ceiling %.3f cyc/%s -- %s\n", ceiling_.cycles_per_unit,
                  unit_name_.c_str(), ceiling_.rationale.c_str());

    const bool have = !results_.empty() && results_[0].counters.available;
    std::printf("\n  %-34s %9s %6s", "variant", "ns/unit", "IQR");
    if (have)
      std::printf(" %9s %6s %9s %9s", "cyc/unit", "IPC", "brmiss/u",
                  "L1miss/u");
    if (!ceiling_.empty() && have)
      std::printf(" %8s", "%ceil");
    std::printf("\n");

    for (const Result &r : results_) {
      std::printf("  %-34s %9.2f %5.1f%%", r.name.c_str(),
                  r.min_s * 1e9 / units_, 100.0 * r.iqr_frac);
      if (have) {
        const double cyc = r.counters.per_unit(units_);
        std::printf(" %9.2f %6.2f %9.4f %9.4f", cyc, r.counters.ipc(),
                    static_cast<double>(r.counters.branch_misses) / units_,
                    static_cast<double>(r.counters.l1d_load_misses) / units_);
        if (!ceiling_.empty())
          std::printf(" %7.1f%%", 100.0 * ceiling_.achieved(cyc));
      }
      std::printf("\n");
    }
    if (!have && !results_.empty())
      std::printf("  (hardware counters unavailable: %s)\n",
                  PerfCounters().reason().c_str());
  }

  std::string title_;
  double units_ = 1;
  std::string unit_name_;
  Ceiling ceiling_;
  std::vector<std::string> notes_;
  std::vector<Variant> variants_;
  std::vector<Result> results_;
  double ghz_ = 0;
  bool ghz_exact_ = false;
  double drift_ = 0;
  std::uint64_t sink_ = 0;
};

} // namespace pte::bench
