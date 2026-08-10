// The device champion meet, against a host oracle written here.
//
// This file is a .cc, which is the mechanical proof of the seam: a test that
// cannot contain device code cannot reach past champion.h, and it cannot link
// pte_cpu either. The oracle below is therefore a *third* implementation of the
// same question -- naive O(D^2) with a std::map of exact keys -- sharing no
// collision-detection code with either the kernel or cpu/champion's reference.
// That is the arrangement REVIEW.md D5 protects: the value of a cross-check is
// that the two sides cannot agree by sharing a bug.
//
// The fixtures the plan names (drafts/gpu_implementation_plan.md 6, "for S-E
// there is no oracle either"):
//
//   * a planted collision at champion widths, injected and required back;
//   * agreement with a brute-force host meet on inputs small enough to have
//     one;
//   * plus, because this path *is* a filter: the exact re-check must reject
//     every false positive, and no window may overflow silently.
#include "gpu/champion/champion.h"

#include "core/uint.h"
#include "gpu/core/device.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

namespace pte::gpu {
namespace {

// The oracle: every pair, an exact key, a map. O(D^2) in time and memory, so it
// is usable only at test sizes -- which is the point of it.
std::vector<core::PairCollision> HostOracle(const std::vector<core::U128> &nu) {
  struct Less {
    bool operator()(const core::U256 &a, const core::U256 &b) const {
      return core::cmp(a, b) < 0;
    }
  };
  std::map<core::U256, std::vector<std::pair<std::uint32_t, std::uint32_t>>,
           Less>
      by_key;
  const std::uint32_t d = static_cast<std::uint32_t>(nu.size());
  std::vector<core::U256> sq(d);
  for (std::uint32_t i = 0; i < d; ++i)
    sq[i] = core::square(nu[i]);
  for (std::uint32_t i = 0; i < d; ++i)
    for (std::uint32_t j = i + 1; j < d; ++j)
      by_key[core::add(sq[i], sq[j])].emplace_back(i, j);

  std::vector<core::PairCollision> out;
  for (const auto &[key, pairs] : by_key) {
    (void)key;
    for (std::size_t p = 0; p + 1 < pairs.size(); ++p)
      for (std::size_t q = p + 1; q < pairs.size(); ++q) {
        const auto &u = pairs[p];
        const auto &v = pairs[q];
        if (u.first == v.first || u.first == v.second || u.second == v.first ||
            u.second == v.second)
          continue;
        out.push_back(core::PairCollision{u.first, u.second, v.first, v.second}
                          .canonical());
      }
  }
  std::sort(out.begin(), out.end());
  return out;
}

// 65 = 1^2 + 8^2 = 4^2 + 7^2, shifted to champion width: nu ~ 2^122 and the key
// ~2^244, so a key compared in 64 or 128 bits would not find it.
struct Planted {
  std::vector<core::U128> nu;
  std::array<std::uint32_t, 4> where;
};

Planted MakePlanted(std::size_t d, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<core::U128> nu;
  std::set<core::U128> used;
  while (nu.size() < d) {
    const core::U128 v = core::u128_from(rng() | 1ull, rng() >> 6);
    if (used.insert(v).second)
      nu.push_back(v);
  }
  const std::uint64_t base[4] = {1, 8, 4, 7};
  std::set<std::uint32_t> slots;
  while (slots.size() < 4)
    slots.insert(static_cast<std::uint32_t>(rng() % d));
  std::array<std::uint32_t, 4> where{};
  int k = 0;
  for (const std::uint32_t slot : slots) {
    where[static_cast<std::size_t>(k)] = slot;
    nu[slot] = core::shl(core::u128_from(base[k]), 119);
    ++k;
  }
  return Planted{nu, where};
}

class ChampionGpuTest : public ::testing::Test {
protected:
  void SetUp() override {
    if (device_count() == 0)
      GTEST_SKIP() << "no CUDA device";
  }

  // Small tables, so the tests exercise many windows and a real load factor
  // rather than one launch into a 512 MB table.
  static ChampionParams SmallTable(int table_log = 12) {
    ChampionParams p;
    p.table_log = table_log;
    p.max_load = 0.5;
    return p;
  }
};

TEST_F(ChampionGpuTest, AgreesWithTheHostOracle) {
  std::mt19937_64 rng(555);
  std::vector<core::U128> nu;
  std::set<core::U128> used;
  // Narrow values on purpose: small keys collide often, so the oracle has
  // something to compare beyond "both found nothing".
  while (nu.size() < 512) {
    const core::U128 v = core::u128_from((rng() % 4096) * 2 + 1);
    if (used.insert(v).second)
      nu.push_back(v);
  }
  ChampionStats stats;
  const auto got = champion_meet(std::span<const core::U128>(nu), SmallTable(),
                                 &stats);
  const auto want = HostOracle(nu);
  EXPECT_FALSE(want.empty()) << "the fixture must contain collisions";
  EXPECT_EQ(got, want);
  EXPECT_EQ(stats.pairs, 512ull * 511 / 2) << "pairs were lost or duplicated";
  EXPECT_EQ(stats.hits, got.size());
}

TEST_F(ChampionGpuTest, RecoversAPlantedCollisionAtChampionWidths) {
  const Planted p = MakePlanted(4096, 20260810);
  ChampionStats stats;
  const auto got =
      champion_meet(std::span<const core::U128>(p.nu), SmallTable(), &stats);

  const core::PairCollision want =
      core::PairCollision{p.where[0], p.where[1], p.where[2], p.where[3]}
          .canonical();
  EXPECT_NE(std::find(got.begin(), got.end(), want), got.end())
      << "the planted collision was not recovered";
  EXPECT_EQ(stats.pairs, 4096ull * 4095 / 2);
  EXPECT_GT(stats.windows, 100u) << "the small table should force many windows";

  // Everything returned is a real collision: the host re-check is the arbiter,
  // and this asserts it did its job rather than passing candidates through.
  for (const core::PairCollision &h : got) {
    EXPECT_TRUE(core::add(core::square(p.nu[h.a]), core::square(p.nu[h.b])) ==
                core::add(core::square(p.nu[h.c]), core::square(p.nu[h.d])));
  }
}

// The answer must not depend on how the key range was cut, on how big the table
// was, or on the nondeterministic order the device appends in. A pair lost at a
// window boundary would show up here and nowhere else.
TEST_F(ChampionGpuTest, AnswerIsIndependentOfTableSize) {
  const Planted p = MakePlanted(1024, 99);
  const auto reference =
      champion_meet(std::span<const core::U128>(p.nu), SmallTable(16));
  ASSERT_FALSE(reference.empty());
  for (const int table_log : {10, 12, 14, 18}) {
    ChampionStats stats;
    const auto got = champion_meet(std::span<const core::U128>(p.nu),
                                   SmallTable(table_log), &stats);
    EXPECT_EQ(got, reference) << "table_log " << table_log
                              << " changed the answer";
    EXPECT_EQ(stats.pairs, 1024ull * 1023 / 2);
  }
}

TEST_F(ChampionGpuTest, RepeatedRunsAgreeExactly) {
  const Planted p = MakePlanted(1024, 4);
  const auto first = champion_meet(std::span<const core::U128>(p.nu),
                                   SmallTable());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(champion_meet(std::span<const core::U128>(p.nu), SmallTable()),
              first)
        << "atomic append is unordered, so the output must be sorted before it "
           "is returned";
}

// A group of g equal keys owes C(g,2) reports, and the "exactly one of two
// equal keys sees the other" argument in champion_impl.cu is what delivers
// that. Three points sharing a key is the smallest case where a naive "report
// the first match and stop" kernel gives the wrong count.
TEST_F(ChampionGpuTest, ReportsEveryPairOfAnEqualKeyGroup) {
  // 50 = 1 + 49 = 25 + 25 (excluded, not disjoint) ; use 325 = 1^2+18^2 =
  // 6^2+17^2 = 10^2+15^2, three representations, so three pairs share a key
  // and C(3,2) = 3 collisions are owed.
  const std::uint64_t v[6] = {1, 18, 6, 17, 10, 15};
  std::vector<core::U128> nu;
  for (const std::uint64_t x : v)
    nu.push_back(core::u128_from(x));
  const auto got = champion_meet(std::span<const core::U128>(nu), SmallTable());
  EXPECT_EQ(got, HostOracle(nu));
  EXPECT_EQ(got.size(), 3u) << "three representations owe C(3,2) = 3 hits";
}

// ---------------------------------------------------------------------------
// Contract: the failures that must be loud
// ---------------------------------------------------------------------------

TEST_F(ChampionGpuTest, RejectsRepeatedNuValues) {
  std::vector<core::U128> nu{core::u128_from(3), core::u128_from(5),
                            core::u128_from(7), core::u128_from(5)};
  EXPECT_THROW(champion_meet(std::span<const core::U128>(nu)),
               std::invalid_argument);
}

TEST_F(ChampionGpuTest, RejectsAMalformedTableSize) {
  const Planted p = MakePlanted(64, 1);
  ChampionParams bad = SmallTable();
  bad.table_log = 0;
  EXPECT_THROW(champion_meet(std::span<const core::U128>(p.nu), bad),
               std::invalid_argument);
  bad = SmallTable();
  bad.max_load = 1.5;
  EXPECT_THROW(champion_meet(std::span<const core::U128>(p.nu), bad),
               std::invalid_argument);
}

// Truncating the candidate list would drop solutions with no diagnostic, so it
// throws instead. This is the same contract DeviceAppendBuffer states: the
// counter keeps rising past capacity precisely so the caller can tell.
TEST_F(ChampionGpuTest, ThrowsRatherThanTruncateCandidates) {
  std::mt19937_64 rng(31337);
  std::vector<core::U128> nu;
  std::set<core::U128> used;
  while (nu.size() < 400) {
    const core::U128 v = core::u128_from((rng() % 800) * 2 + 1);
    if (used.insert(v).second)
      nu.push_back(v);
  }
  ChampionParams tiny = SmallTable();
  tiny.max_candidates = 1;
  EXPECT_THROW(champion_meet(std::span<const core::U128>(nu), tiny),
               std::runtime_error);
}

TEST_F(ChampionGpuTest, FewerThanFourPointsHasNoSolutions) {
  std::vector<core::U128> nu{core::u128_from(3), core::u128_from(5),
                            core::u128_from(7)};
  EXPECT_TRUE(champion_meet(std::span<const core::U128>(nu)).empty());
}

} // namespace
} // namespace pte::gpu
