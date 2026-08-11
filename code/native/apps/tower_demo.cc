// G1: the local model, measured.
//
//   ./tower_demo [max_modulus] [max_depth] [only_p]
//
// Prints three tables, in the order drafts/local_model_plan.md needs them:
//
//   1. V_p(k,d), the least node-valuation cap admitting a depth-d tower over
//      Z/p^k. This is the quantity the whole plan turns on. V = k means
//      "nothing nondegenerate is visible at this precision" -- which grows into
//      a nonexistence proof if it keeps growing with k, and means *nothing at
//      all* below the precision a real solution needs.
//
//   2. The exemplar control. The section-7 E_4 is a real depth-4 tower, so its
//      reduction must appear once k is large enough; the cap at which it does
//      is ht_p, its largest node valuation over Z. Without this column a null
//      in table 1 is indistinguishable from insufficient precision.
//
//   3. The translation-orbit structure, which is what decides how far this can
//      be pushed: if each level's state set is a union of complete orbits, a
//      larger run stores O(1) representatives rather than O(p^k) states.
//
// The default budget reaches 3^7 = 2187, which is not an arbitrary stopping
// point: paper.md 4's C'_16 = 2^11 3^6 5^4 7^3 ... forces k >= 7 at p = 3 for
// any depth-5 claim to mean anything, and this is the first precision that
// clears it.
#include "cpu/tower/tower.h"
#include "cpu/tower/tower_oracle.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using pte::cpu::TowerCounts;
using pte::cpu::TowerOracle;
using pte::cpu::TowerParams;
using pte::cpu::TowerTooLarge;
using pte::cpu::ValuationFloor;

constexpr std::uint32_t kAllPrimes[] = {2, 3, 5, 7, 11, 13};

// Set from argv[3] to sweep one prime deeply rather than all of them shallowly,
// which is what pushing past paper.md 4's precision floor requires.
std::vector<std::uint32_t> kPrimes(std::begin(kAllPrimes), std::end(kAllPrimes));

double Now() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

// p^k, or 0 on overflow past `budget`.
std::uint64_t Pow(std::uint32_t p, std::uint32_t k, std::uint64_t budget) {
  std::uint64_t m = 1;
  for (std::uint32_t i = 0; i < k; ++i) {
    m *= p;
    if (m > budget)
      return 0;
  }
  return m;
}

std::string Join(const std::vector<std::uint64_t> &values) {
  std::string out;
  for (const std::uint64_t v : values) {
    if (!out.empty())
      out += ' ';
    out += std::to_string(v);
  }
  return out;
}

// V_p(k,d) only ever needs existence, and turning counting off is what buys the
// shortcut in valuation_floor: one run at cap = k-1 settles the common "nothing
// at this precision" case instead of k ascending runs.
ValuationFloor Floor(std::uint32_t p, std::uint32_t k, std::uint32_t depth) {
  TowerParams base;
  base.count_towers = false;
  return pte::cpu::valuation_floor(p, k, depth, base);
}

void PrintFloorTable(std::uint64_t budget, std::uint32_t max_depth) {
  std::printf(
      "\n1. V_p(k,d): least node-valuation cap admitting a depth-d tower over "
      "Z/p^k\n");
  std::printf("   (V = k, shown as k*, means nothing nondegenerate is visible "
              "at that precision)\n\n");
  std::printf("   %-12s", "p^k");
  for (std::uint32_t d = 2; d <= max_depth; ++d)
    std::printf("  d=%-4u", d);
  std::printf("   %-28s %8s\n", "states/level (uncapped, max d)", "seconds");

  for (const std::uint32_t p : kPrimes) {
    bool any = false;
    for (std::uint32_t k = 1;; ++k) {
      const std::uint64_t m = Pow(p, k, budget);
      if (m == 0)
        break;
      any = true;
      const double t0 = Now();
      const std::string label =
          std::to_string(p) + "^" + std::to_string(k) + "=" + std::to_string(m);
      std::printf("   %-12s", label.c_str());
      std::fflush(stdout);
      // A row that does not fit is a measurement about the cost model rather
      // than a crash: the state table is what bounds this, and *which* (p,k,d)
      // hits the wall is the datum drafts/local_model_plan.md 5.1 argues from.
      bool oversized = false;
      for (std::uint32_t d = 2; d <= max_depth && !oversized; ++d) {
        try {
          const ValuationFloor v = Floor(p, k, d);
          const std::string cell =
              std::to_string(v.floor) + (v.nondegenerate() ? "" : "*");
          std::printf("  %-6s", cell.c_str());
        } catch (const TowerTooLarge &) {
          oversized = true;
        }
        std::fflush(stdout);
      }
      if (oversized) {
        std::printf("   over max_states -- see plan 5.1%14.2f\n", Now() - t0);
        continue;
      }
      TowerParams params;
      params.p = p;
      params.k = k;
      params.depth = max_depth;
      params.count_towers = false;
      const TowerCounts counts = pte::cpu::tower_counts(params);
      std::printf("   %-28s %8.2f\n", Join(counts.states).c_str(), Now() - t0);
    }
    if (any)
      std::printf("\n");
  }
}

void PrintExemplarTable(std::uint64_t budget) {
  std::printf(
      "2. The control: the section-7 exemplar E_4, reduced mod p^k\n"
      "   Its least admitting cap is min(ht_p, k), so the column stops rising\n"
      "   exactly at ht_p -- recovered from the reduced tuple alone, with no\n"
      "   big-integer arithmetic. Table 1 cannot be read below these numbers.\n\n");
  for (const std::uint32_t p : kPrimes) {
    std::vector<std::uint32_t> floors;
    std::uint32_t last_k = 0;
    for (std::uint32_t k = 1;; ++k) {
      const std::uint64_t m = Pow(p, k, budget);
      if (m == 0)
        break;
      const TowerOracle oracle(p, k);
      floors.push_back(
          oracle.tuple_valuation_floor(pte::cpu::exemplar_e4_mod(oracle.modulus())));
      last_k = k;
    }
    if (floors.empty())
      continue;
    std::string row;
    for (const std::uint32_t f : floors)
      row += std::to_string(f) + " ";
    const std::string verdict =
        floors.back() < last_k
            ? "-> ht_" + std::to_string(p) + " = " + std::to_string(floors.back())
            : "-> ht_p is still >= k here; raise the budget";
    std::printf("   p=%-3u k=1..%-2u  %-26s %s\n", p, last_k, row.c_str(),
                verdict.c_str());
  }
  std::printf("\n");
}

void PrintOrbitTable(std::uint64_t budget, std::uint32_t max_depth) {
  std::printf(
      "3. Translation orbits. Every state at level >= 1 has the form S - L, so\n"
      "   each level's state set should be a union of complete orbits: misses\n"
      "   counts rotations that are not themselves states, and must be 0.\n"
      "   Where it holds, states ~ orbits * p^k (exactly, minus the stabilised\n"
      "   orbits, plus the absorbing empty state) and the object worth storing\n"
      "   is the representative, not the state.\n\n");
  std::printf("   %-12s %-30s %-16s %s\n", "p^k", "states/level", "orbits/level",
              "misses");
  for (const std::uint32_t p : kPrimes)
    for (std::uint32_t k = 1;; ++k) {
      // The diagnostic costs an extra states * m * m/64 pass, so it runs on a
      // tighter budget than the table it explains.
      const std::uint64_t m = Pow(p, k, budget);
      if (m == 0)
        break;
      TowerParams params;
      params.p = p;
      params.k = k;
      params.depth = max_depth;
      params.count_towers = false;
      params.count_orbits = true;
      TowerCounts counts;
      try {
        counts = pte::cpu::tower_counts(params);
      } catch (const TowerTooLarge &) {
        continue;
      }
      std::uint64_t misses = 0;
      for (const std::uint64_t x : counts.orbit_misses)
        misses += x;
      const std::string label =
          std::to_string(p) + "^" + std::to_string(k) + "=" + std::to_string(m);
      std::printf("   %-12s %-30s %-16s %llu\n", label.c_str(),
                  Join(counts.states).c_str(), Join(counts.orbits).c_str(),
                  static_cast<unsigned long long>(misses));
    }
  std::printf("\n");
}

} // namespace

int main(int argc, char **argv) {
  const std::uint64_t budget = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 2400;
  const std::uint32_t max_depth =
      argc > 2 ? static_cast<std::uint32_t>(std::strtoul(argv[2], nullptr, 10)) : 5;
  if (argc > 3)
    kPrimes = {static_cast<std::uint32_t>(std::strtoul(argv[3], nullptr, 10))};

  std::printf("G1 local model -- modulus budget p^k <= %llu, depth <= %u\n",
              static_cast<unsigned long long>(budget), max_depth);
  std::printf("paper.md 4's C'_16 forces the precision at which a depth-5 claim\n"
              "first means anything: k >= 12, 7, 5, 4, 3, 3, 2 at p = 2, 3, 5,\n"
              "7, 11, 13, <= 53. Below that a null is not a result.\n");

  PrintFloorTable(budget, max_depth);
  PrintExemplarTable(budget);
  PrintOrbitTable(budget < 400 ? budget : 400, max_depth);
  return 0;
}
