// The property the whole of S-E rests on: value-windowing enumerates every pair
// exactly once.
//
// This is check 4 of drafts/verify_gpu_plan.py, moved into the build. The
// Python version proved the idea; this one guards the implementation, and it is
// the more important of the two because the failure it looks for is invisible
// in output:
//
//   * a pair enumerated **twice** is wasted work and duplicate hits --
//   annoying,
//     and immediately visible;
//   * a pair enumerated **zero** times is a lost solution, reported as "no
//     solution here", forever.
//
// So every test below checks coverage and uniqueness together, against a
// brute-force set of all C(D,2) pairs.
#include "core/window.h"

#include "core/uint.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace pte::core {
namespace {

// A sorted square array standing in for a point set's sq[]. Values are spread
// over the full 244-bit champion range when `wide` is set, and kept narrow
// otherwise, because the two exercise different limbs of every comparison.
std::vector<U256> MakeSq(std::size_t d, std::uint64_t seed, bool wide) {
  std::mt19937_64 rng(seed);
  std::vector<U128> nu;
  nu.reserve(d);
  for (std::size_t i = 0; i < d; ++i)
    nu.push_back(wide ? u128_from(rng(), rng() >> 6) : u128_from(rng() >> 40));
  std::vector<U256> sq;
  sq.reserve(d);
  for (const U128 &v : nu)
    sq.push_back(square(v));
  std::sort(sq.begin(), sq.end(),
            [](const U256 &a, const U256 &b) { return cmp(a, b) < 0; });
  return sq;
}

struct Enumerated {
  std::multiset<std::pair<std::uint32_t, std::uint32_t>> pairs;
  std::uint64_t windows = 0;
  std::uint64_t max_window_pairs = 0;
};

// Runs the planner and records every (i, j) it produces, plus a check that each
// key really lies inside the window it was produced in -- the property that
// makes "exactly one window" true in the first place.
Enumerated RunAll(std::span<const U256> sq, std::uint64_t cap) {
  Enumerated out;
  WindowScratch scratch(sq.size());
  WindowRanges ranges = scratch.ranges();
  for_each_window(
      sq, cap, ranges, [&](const KeyWindow &w, const WindowRanges &r) {
        ++out.windows;
        out.max_window_pairs = std::max(out.max_window_pairs, r.pairs);
        std::uint64_t seen = 0;
        for (std::size_t i = 0; i < sq.size(); ++i) {
          for (std::uint32_t j = r.lo[i]; j < r.hi[i]; ++j) {
            EXPECT_GT(j, i) << "a pair must be enumerated with i<j";
            const U256 key = add(sq[i], sq[j]);
            EXPECT_GE(cmp(key, w.lo), 0) << "key below its window";
            EXPECT_LT(cmp(key, w.hi), 0) << "key above its window";
            out.pairs.emplace(static_cast<std::uint32_t>(i), j);
            ++seen;
          }
        }
        EXPECT_EQ(seen, r.pairs) << "pairs count disagrees with the "
                                    "ranges it describes";
      });
  return out;
}

void ExpectExactCover(const Enumerated &got, std::size_t d) {
  std::multiset<std::pair<std::uint32_t, std::uint32_t>> want;
  for (std::uint32_t i = 0; i < d; ++i)
    for (std::uint32_t j = i + 1; j < d; ++j)
      want.emplace(i, j);
  ASSERT_EQ(got.pairs.size(), want.size())
      << "enumerated " << got.pairs.size() << " pairs, want " << want.size();
  EXPECT_TRUE(got.pairs == want) << "the enumerated pair set is not C(D,2)";
}

void ExpectExactCoverAtEveryCap(std::size_t d,
                                std::initializer_list<std::uint64_t> caps) {
  for (const bool wide : {false, true}) {
    const auto sq = MakeSq(d, 1234 + d, wide);
    for (const std::uint64_t cap : caps) {
      const Enumerated got = RunAll(sq, cap);
      ExpectExactCover(got, d);
      EXPECT_LE(got.max_window_pairs, cap)
          << "a window exceeded the cap (D=" << d << ", cap=" << cap << ")";
    }
  }
}

// A cap of 1 is the extreme case -- one window per distinct key value, so the
// planner bisects all the way down and every leaf is a single pair. It is the
// cheapest way to exercise the recursion to its floor, so it runs at the small
// sizes; the larger D only needs caps that leave the tree shallow.
TEST(Window, EnumeratesEveryPairExactlyOnce) {
  for (const std::size_t d :
       {std::size_t{2}, std::size_t{3}, std::size_t{7}, std::size_t{64}})
    ExpectExactCoverAtEveryCap(d, {1, 3, 17, 1000, 1000000});
  ExpectExactCoverAtEveryCap(400, {17, 1000, 1000000});
}

// The cap is the caller's memory bound: it is what sizes the device hash table,
// so exceeding it silently would overflow that table rather than merely run
// slowly. A tighter cap must only ever mean more windows.
TEST(Window, TighterCapsMeanMoreWindowsAndNeverMorePairs) {
  const auto sq = MakeSq(300, 99, true);
  std::uint64_t prev_windows = 0;
  for (const std::uint64_t cap : {100000u, 10000u, 1000u, 100u, 10u}) {
    const Enumerated got = RunAll(sq, cap);
    ExpectExactCover(got, 300);
    EXPECT_GE(got.windows, prev_windows);
    prev_windows = got.windows;
  }
}

// A duplicated sq value is out of contract for the meet (nu is injective on a
// canonized point set), but the planner must still fail loudly rather than
// hand back a window it cannot cut.
TEST(Window, ThrowsWhenOneKeyValueExceedsTheCap) {
  std::vector<U256> sq(100, square(u128_from(12345)));
  WindowScratch scratch(sq.size());
  WindowRanges ranges = scratch.ranges();
  EXPECT_THROW(for_each_window(std::span<const U256>(sq), 10, ranges,
                               [](const KeyWindow &, const WindowRanges &) {}),
               WindowIndivisible);
}

TEST(Window, RejectsUnsortedInput) {
  std::vector<U256> sq{square(u128_from(9)), square(u128_from(4))};
  WindowScratch scratch(sq.size());
  WindowRanges ranges = scratch.ranges();
  EXPECT_THROW(for_each_window(std::span<const U256>(sq), 100, ranges,
                               [](const KeyWindow &, const WindowRanges &) {}),
               std::invalid_argument);
}

TEST(Window, FewerThanTwoPointsProducesNoWindows) {
  for (const std::size_t d : {std::size_t{0}, std::size_t{1}}) {
    const auto sq = MakeSq(d, 5, false);
    WindowScratch scratch(std::max<std::size_t>(d, 1));
    WindowRanges ranges = scratch.ranges();
    int calls = 0;
    for_each_window(std::span<const U256>(sq), 100, ranges,
                    [&](const KeyWindow &, const WindowRanges &) { ++calls; });
    EXPECT_EQ(calls, 0);
  }
}

// The two-pointer sweep is an optimization of "binary search the j-range";
// this checks it against the definition directly, since a sweep that loses
// monotonicity drops pairs on one side of a boundary only.
TEST(Window, RangesMatchTheDefinition) {
  const auto sq = MakeSq(200, 31337, true);
  const U256 lo = add(sq[10], sq[40]);
  const U256 hi = add(sq[60], sq[190]);
  WindowScratch scratch(sq.size());
  WindowRanges ranges = scratch.ranges();
  window_ranges(std::span<const U256>(sq), KeyWindow{lo, hi}, ranges);

  std::uint64_t want_pairs = 0;
  for (std::size_t i = 0; i < sq.size(); ++i) {
    std::size_t want_lo = sq.size(), want_hi = sq.size();
    for (std::size_t j = i + 1; j < sq.size(); ++j) {
      const U256 key = add(sq[i], sq[j]);
      if (cmp(key, lo) >= 0 && want_lo == sq.size())
        want_lo = j;
      if (cmp(key, hi) >= 0 && want_hi == sq.size()) {
        want_hi = j;
        break;
      }
    }
    if (want_lo > want_hi)
      want_lo = want_hi;
    EXPECT_EQ(ranges.lo[i], want_lo) << "lo range at i=" << i;
    EXPECT_EQ(ranges.hi[i], want_hi) << "hi range at i=" << i;
    want_pairs += want_hi - want_lo;
  }
  EXPECT_EQ(ranges.pairs, want_pairs);
}

} // namespace
} // namespace pte::core
