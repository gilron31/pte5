// The machine model a benchmark is measured against.
//
// A benchmark that reports only "42 ns/norm" tells you nothing about whether to
// keep working on it. The number that decides is *percent of ceiling*: 90% means
// stop, 10% means there is an order of magnitude on the table. So every bench in
// this directory declares the theoretical cost of the work it does, and the
// harness prints achieved against it.
//
// ---------------------------------------------------------------------------
// The machine (AMD Ryzen 7 8845HS, Zen 4)
// ---------------------------------------------------------------------------
//
// Only the parameters the enumeration actually leans on. Sources are AMD's
// Software Optimization Guide for Zen 4 (PUB 57647) and Agner Fog's tables;
// where they disagree the more conservative figure is used, because a ceiling
// that is too generous makes real code look worse than it is.
//
//   integer ALU pipes            4    (ALU0-3; simple add/sub/cmp/logic)
//   load ports                   3    (2x256-bit + 1, but 3 scalar loads/cycle)
//   store ports                  2
//   64-bit imul                  1/cycle throughput, 3 cycle latency (ALU1)
//   branch                       2 taken/cycle predicted, ~13 cycle mispredict
//   L1d                          32 KB, 4-5 cycle load-to-use
//   L2                           1 MB per core, ~14 cycles
//   L3                           16 MB shared (8845HS), ~50 cycles
//
// The frequency is deliberately NOT a constant here. It is measured per run by
// the harness dynamometer, because on this laptop it ranges 0.4-5.1 GHz and a
// hardcoded clock is how a 2x error gets into a table.
#pragma once

#include <string>
#include <utility>

namespace pte::bench {

// Zen 4 issue limits, per core per cycle.
struct Machine {
  static constexpr double kIntAluPipes = 4.0;
  static constexpr double kLoadPorts = 3.0;
  static constexpr double kStorePorts = 2.0;
  static constexpr double kImul64PerCycle = 1.0;
  static constexpr double kBranchMispredictCycles = 13.0;
  static constexpr int kL1dBytes = 32 * 1024;
  static constexpr int kL2Bytes = 1024 * 1024;
  static constexpr int kL3Bytes = 16 * 1024 * 1024;
};

// The floor, in cycles, for one unit of a bench's work -- plus why.
//
// The rationale string is not decoration. A ceiling with no stated derivation
// is unfalsifiable, and these are estimates: they assume perfect scheduling, no
// dependency stalls and no misses, which is exactly why hitting 30-40% of one is
// normal for real code rather than a scandal.
struct Ceiling {
  double cycles_per_unit = 0;
  std::string rationale;

  bool empty() const { return cycles_per_unit <= 0; }

  // What fraction of the ceiling a measurement achieved.
  double achieved(double measured_cycles_per_unit) const {
    return measured_cycles_per_unit > 0
               ? cycles_per_unit / measured_cycles_per_unit
               : 0.0;
  }
};

// A ceiling from a plain op count issued at the ALU width, i.e. the best case
// for straight-line integer work with no memory pressure.
inline Ceiling AluBound(double ops_per_unit, std::string why) {
  return Ceiling{ops_per_unit / Machine::kIntAluPipes, std::move(why)};
}

// A ceiling from load count, for work whose limit is feeding the ALUs.
inline Ceiling LoadBound(double loads_per_unit, std::string why) {
  return Ceiling{loads_per_unit / Machine::kLoadPorts, std::move(why)};
}

// The larger of two constraints -- use when a kernel is plausibly limited by
// either, which is the honest position for most of this code.
inline Ceiling Max(const Ceiling &a, const Ceiling &b) {
  return a.cycles_per_unit >= b.cycles_per_unit ? a : b;
}

} // namespace pte::bench
