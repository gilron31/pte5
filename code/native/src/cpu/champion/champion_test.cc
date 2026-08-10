// S-E has no oracle, which is the whole difficulty of testing it: there is no
// catalogue of champion collisions to compare against, and at D = 2^19 there is
// no brute force either. So the fixtures are the two the plan names
// (drafts/gpu_implementation_plan.md 6), plus the invariant that makes the
// windowing checkable on its own:
//
//   1. **Small-D against the production meet.** On real point sets, where nu is
//      64-bit and D is small, `champion_meet` must agree exactly with
//      `pte::cpu::meet` and `meet_oracle` -- three implementations of one
//      question, sharing no collision-detection code. This is the fixture that
//      would catch a wrong answer in the ordinary regime.
//   2. **A planted collision at champion scale.** Inject a known equal-key
//      quadruple into a large synthetic point set with 122-bit nu values and
//      require it back. This is the only fixture that exercises the widths S-E
//      actually runs at, and the only one that would catch a 256-bit key being
//      compared as 64.
//   3. **Pair conservation.** stats.pairs must equal C(D,2) exactly, on every
//      input and at every window cap. A lost pair is a lost solution and has no
//      other symptom.
#include "cpu/champion/champion.h"

#include "core/window.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/gaussian/gaussian.h"
#include "cpu/meet/meet.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

namespace pte::cpu {
namespace {

std::vector<core::U128> ToU128(const std::vector<norm_t> &nu) {
  std::vector<core::U128> out;
  out.reserve(nu.size());
  for (const norm_t v : nu)
    out.push_back(core::u128_from(static_cast<std::uint64_t>(v)));
  return out;
}

// The production meet's answer, canonized into the shared hit type so the two
// searches are comparable at all.
std::vector<core::PairCollision> MeetAnswer(const std::vector<norm_t> &nu) {
  std::vector<core::PairCollision> out;
  for (const E4Sol &s :
       meet_collect<exact_key_t>(std::span<const norm_t>(nu))) {
    const E4Sol c = s.canonical();
    out.push_back(core::PairCollision{
        static_cast<std::uint32_t>(c.a), static_cast<std::uint32_t>(c.b),
        static_cast<std::uint32_t>(c.c), static_cast<std::uint32_t>(c.d)}
                      .canonical());
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<core::PairCollision> OracleAnswer(const std::vector<norm_t> &nu) {
  std::vector<core::PairCollision> out;
  for (const E4Sol &s : meet_oracle(nu)) {
    const E4Sol c = s.canonical();
    out.push_back(core::PairCollision{
        static_cast<std::uint32_t>(c.a), static_cast<std::uint32_t>(c.b),
        static_cast<std::uint32_t>(c.c), static_cast<std::uint32_t>(c.d)}
                      .canonical());
  }
  std::sort(out.begin(), out.end());
  return out;
}

// ---------------------------------------------------------------------------
// 1. Small D, real point sets, against the two existing meets
// ---------------------------------------------------------------------------

TEST(ChampionMeet, AgreesWithTheProductionMeetOnTheExemplarNorm) {
  // paper.md's section 7 exemplar: 13 * 89 * 173 * 233 * 709, the norm whose
  // solution the whole project is calibrated on.
  const std::vector<Gauss> points = point_set({13, 89, 173, 233, 709});
  const std::vector<norm_t> nu = nu_values(points);
  ASSERT_EQ(nu.size(), 16u);

  ChampionStats stats;
  const auto got =
      champion_meet(std::span<const core::U128>(ToU128(nu)), {}, &stats);
  const auto want = MeetAnswer(nu);
  EXPECT_FALSE(want.empty()) << "the exemplar must contain a solution";
  EXPECT_EQ(got, want);
  EXPECT_EQ(got, OracleAnswer(nu));
  EXPECT_EQ(stats.pairs, nu.size() * (nu.size() - 1) / 2);
}

TEST(ChampionMeet, AgreesWithTheProductionMeetOverASweepOfNorms) {
  Enumerator enumerator(200000);
  int checked = 0, with_hits = 0;
  enumerator.for_each_norm([&](const DecomposedNorm &norm) {
    const std::vector<norm_t> nu = nu_values(enumerator.point_set(norm));
    if (nu.size() < 4)
      return;
    ++checked;
    const auto got = champion_meet(std::span<const core::U128>(ToU128(nu)));
    const auto want = MeetAnswer(nu);
    if (!want.empty())
      ++with_hits;
    ASSERT_EQ(got, want) << "disagreement at L1 = " << norm.l1
                         << " (D = " << nu.size() << ")";
  });
  EXPECT_GT(checked, 1000) << "the sweep must actually cover norms";
  EXPECT_GT(with_hits, 0) << "a sweep with no solutions proves nothing about "
                             "agreement on solutions";
}

// A tighter cap forces bisection, which is where a pair could be lost at a
// window boundary -- the failure mode with no visible symptom. Same input, many
// partitions, one answer.
TEST(ChampionMeet, AnswerIsIndependentOfTheWindowCap) {
  const std::vector<norm_t> nu = nu_values(point_set({13, 89, 173, 233, 709}));
  const auto reference = champion_meet(std::span<const core::U128>(ToU128(nu)));
  ASSERT_FALSE(reference.empty());

  // Note the floor: a solution *is* two pairs sharing one key value, so a cap
  // of 1 cannot be satisfied on a norm that has one. Bisection cannot separate
  // equal keys, and WindowIndivisible says so rather than dropping a pair --
  // which is the whole reason that error exists.
  {
    ChampionParams too_tight;
    too_tight.max_pairs_per_window = 1;
    EXPECT_THROW(
        champion_meet(std::span<const core::U128>(ToU128(nu)), too_tight),
        core::WindowIndivisible);
  }

  for (const std::uint64_t cap : {2u, 5u, 17u, 120u, 1000u}) {
    ChampionStats stats;
    ChampionParams params;
    params.max_pairs_per_window = cap;
    const auto got =
        champion_meet(std::span<const core::U128>(ToU128(nu)), params, &stats);
    EXPECT_EQ(got, reference) << "cap " << cap << " changed the answer";
    EXPECT_EQ(stats.pairs, nu.size() * (nu.size() - 1) / 2)
        << "cap " << cap << " lost or duplicated pairs";
  }
}

// ---------------------------------------------------------------------------
// 2. Champion scale: the planted collision
// ---------------------------------------------------------------------------

// Build a point set of `d` distinct nu values at champion width, with exactly
// one planted equal-key quadruple.
//
// The plant is 65 = 1^2 + 8^2 = 4^2 + 7^2 -- two representations of the same
// sum of two squares, which is precisely an E_4 configuration -- with each
// value shifted left by 119 bits. That puts nu at ~122 bits and the key at
// 65 * 2^238, i.e. ~244: the widths S-E actually runs at, and past every 32-,
// 64- and 128-bit boundary on the way, so a key compared in any narrower type
// would fail to find it.
struct Planted {
  std::vector<core::U128> nu;
  std::array<std::uint32_t, 4> where; // indices of the planted quadruple
};

Planted MakePlanted(std::size_t d, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<core::U128> nu;
  std::set<core::U128> used;
  // Random 122-bit *odd* values, like a real nu: L1 is odd, so nu = m^2 - n^2
  // is odd. That also keeps every random value distinct from the planted ones,
  // which are even by construction.
  while (nu.size() < d) {
    const core::U128 v = core::u128_from(rng() | 1ull, rng() >> 6);
    if (used.insert(v).second)
      nu.push_back(v);
  }

  const std::uint64_t base[4] = {1, 8, 4, 7};
  std::array<std::uint32_t, 4> where{};
  std::set<std::uint32_t> slots;
  while (slots.size() < 4)
    slots.insert(static_cast<std::uint32_t>(rng() % d));
  int k = 0;
  for (const std::uint32_t slot : slots) {
    where[static_cast<std::size_t>(k)] = slot;
    nu[slot] = core::shl(core::u128_from(base[k]), 119);
    ++k;
  }
  return Planted{nu, where};
}

TEST(ChampionMeet, RecoversAPlantedCollisionAtChampionWidths) {
  const Planted p = MakePlanted(4096, 20260810);
  ChampionStats stats;
  ChampionParams params;
  params.max_pairs_per_window = 1u << 14; // force many windows
  const auto got =
      champion_meet(std::span<const core::U128>(p.nu), params, &stats);

  core::PairCollision want{p.where[0], p.where[1], p.where[2], p.where[3]};
  want = want.canonical();
  EXPECT_NE(std::find(got.begin(), got.end(), want), got.end())
      << "the planted collision was not recovered";
  EXPECT_EQ(stats.pairs, 4096ull * 4095 / 2) << "pairs were lost";
  EXPECT_GT(stats.windows, 100u) << "the cap should have forced many windows";

  // Every reported hit must be a real one: recompute the key exactly.
  for (const core::PairCollision &h : got) {
    const core::U256 left =
        core::add(core::square(p.nu[h.a]), core::square(p.nu[h.b]));
    const core::U256 right =
        core::add(core::square(p.nu[h.c]), core::square(p.nu[h.d]));
    EXPECT_TRUE(left == right) << "a reported hit is not a collision";
  }
}

TEST(ChampionMeet, HitsAreUniqueAndCanonical) {
  const Planted p = MakePlanted(2048, 7);
  const auto got = champion_meet(std::span<const core::U128>(p.nu));
  for (const core::PairCollision &h : got)
    EXPECT_EQ(h, h.canonical());
  EXPECT_TRUE(std::is_sorted(got.begin(), got.end()));
  EXPECT_EQ(std::adjacent_find(got.begin(), got.end()), got.end())
      << "a duplicate hit means a pair was enumerated in two windows";
}

// ---------------------------------------------------------------------------
// 3. Contract
// ---------------------------------------------------------------------------

TEST(ChampionMeet, RejectsRepeatedNuValues) {
  // nu is injective on a canonized point set; a repeat means the caller passed
  // raw products, and the disjointness argument no longer holds.
  std::vector<core::U128> nu{core::u128_from(3), core::u128_from(5),
                             core::u128_from(7), core::u128_from(5)};
  EXPECT_THROW(champion_meet(std::span<const core::U128>(nu)),
               std::invalid_argument);
}

TEST(ChampionMeet, FewerThanFourPointsHasNoSolutions) {
  std::vector<core::U128> nu{core::u128_from(3), core::u128_from(5),
                             core::u128_from(7)};
  EXPECT_TRUE(champion_meet(std::span<const core::U128>(nu)).empty());
  EXPECT_TRUE(champion_meet(std::span<const core::U128>()).empty());
}

} // namespace
} // namespace pte::cpu
