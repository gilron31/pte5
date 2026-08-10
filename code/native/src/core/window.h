// Value-windowing: the partition of a champion norm's pair space that makes
// S-E affordable.
//
// ---------------------------------------------------------------------------
// The problem, and the 1000x that fixes it
// ---------------------------------------------------------------------------
//
// A champion norm has D ~ 2^19 points, so C(D,2) = 1.374e11 pair keys. At 16
// bytes each that is 2.2 TB: the keys cannot be materialised, and a collision
// search has to process them in slices that fit in memory.
//
// Both earlier drafts sliced by *hash*: keep only the keys whose hash lands in
// a residue window of width 1/R, and repeat R times. A hash window is
// unpredictable, so that recomputes every key R times -- R * C(D,2) = 1.4e14
// key evaluations at R = 10^3.
//
// Slicing by key *value* instead costs each key once. With sq[] sorted
// ascending, sq[i] + sq[j] is increasing in j, so the pairs whose key falls in
// [lo, hi) are, for each i, a **contiguous run of j** -- and every pair belongs
// to exactly one window, because the windows are a half-open partition of an
// interval containing every key:
//
//     total work = C(D,2) + R*D   against   R * C(D,2)
//
// That is drafts/gpu_implementation_plan.md 3.1, and it is three orders of
// magnitude for a change of *which* quantity the window is cut on.
//
// ---------------------------------------------------------------------------
// Why the comparisons are exact, and what happens if they are not
// ---------------------------------------------------------------------------
//
// A key is ~244 bits (core/u256.h). If the comparison against a window bound
// were approximate -- a double, or a truncation to 64 bits -- a pair could be
// placed in a window on one side of the boundary while its equal-key partner
// lands on the other, and the two are then never compared. **No output changes,
// no diagnostic appears, and the solution is gone.** So every comparison here
// is a full 256-bit integer compare, and this file has no floating point in it
// at all.
//
// Truncation *is* sound in the other direction, inside a window: once two pairs
// are in the same window, comparing them mod 2^64 has false positives but no
// false negatives, which is the same asymmetry cpu/core/precision.h relies on.
// The device kernel uses that; the windowing may not.
//
// ---------------------------------------------------------------------------
// Sizing the windows
// ---------------------------------------------------------------------------
//
// The caller states a cap on pairs per window (whatever its collision table
// holds) rather than a window count, and the planner bisects until every leaf
// fits. Bisection keeps the exact-partition property for free: [lo,mid) and
// [mid,hi) are disjoint and cover [lo,hi). It also adapts to a lopsided key
// distribution, which a fixed R does not -- and champion keys are lopsided,
// since sq[] spans a wide range and the pair sums pile up near the top.
#pragma once

#include "core/uint.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace pte::core {

// A half-open interval of key values, [lo, hi).
struct KeyWindow {
  U256 lo, hi;
};

// For one window: the j-range [lo[i], hi[i]) of every i, and the total pair
// count.
//
// The ranges are *views into caller-owned storage*, not vectors of their own.
// That is what lets the device path have the planner write straight into pinned
// staging memory instead of filling a std::vector and copying it -- the "the
// producer must write into the pinned buffer, not copy into it" rule from
// REVIEW.md A2(a). See WindowScratch for the plain-memory case.
struct WindowRanges {
  std::span<std::uint32_t> lo, hi;
  std::uint64_t pairs = 0;
};

// Ordinary heap storage for the ranges, for callers with no pinned buffer.
class WindowScratch {
public:
  explicit WindowScratch(std::size_t d) : lo_(d), hi_(d) {}
  WindowRanges ranges() { return WindowRanges{lo_, hi_, 0}; }

private:
  std::vector<std::uint32_t> lo_, hi_;
};

// Fills `out` with the j-ranges of `window` over `sq`, which must be sorted
// ascending. `out.lo` and `out.hi` must each have sq.size() entries.
//
// Two pointers, not two binary searches: the thresholds are non-increasing in i
// (sq[i] grows, so the j that first reaches a fixed bound shrinks), so one
// descending sweep costs O(D) per window against the O(D log D) the plan
// budgets. The `max(i+1, ...)` is what keeps each unordered pair to a single
// (i, j) with i < j.
inline void window_ranges(std::span<const U256> sq, const KeyWindow &window,
                          WindowRanges &out) {
  const std::size_t d = sq.size();
  if (out.lo.size() < d || out.hi.size() < d)
    throw std::invalid_argument("window_ranges: range storage too small");

  std::uint64_t pairs = 0;
  std::size_t t_lo = d, t_hi = d;
  for (std::size_t i = 0; i < d; ++i) {
    const U256 &si = sq[i];
    while (t_lo > 0 && cmp(add(si, sq[t_lo - 1]), window.lo) >= 0)
      --t_lo;
    while (t_hi > 0 && cmp(add(si, sq[t_hi - 1]), window.hi) >= 0)
      --t_hi;
    const std::size_t begin = std::max(i + 1, t_lo);
    const std::size_t end = std::max(begin, t_hi);
    out.lo[i] = static_cast<std::uint32_t>(begin);
    out.hi[i] = static_cast<std::uint32_t>(end);
    pairs += static_cast<std::uint64_t>(end - begin);
  }
  out.pairs = pairs;
}

// Thrown when a single key value carries more pairs than the cap allows, so
// bisection cannot make progress. It is a distinct type because the caller's
// only remedy is a bigger table -- retrying with more windows will not help,
// and silently processing an over-full window would overflow the collision
// table instead.
class WindowIndivisible : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

namespace detail {

template <typename F>
void window_descend(std::span<const U256> sq, const KeyWindow &window,
                    std::uint64_t max_pairs, WindowRanges &scratch, F &f) {
  window_ranges(sq, window, scratch);
  if (scratch.pairs == 0)
    return; // empty windows cost nothing and are not worth a launch
  if (scratch.pairs <= max_pairs) {
    f(window, scratch);
    return;
  }
  // Halving is exact: (hi - lo) >> 1 never rounds away a value, so the two
  // halves still partition the parent.
  const U256 mid = add(window.lo, shr(sub(window.hi, window.lo), 1));
  if (cmp(mid, window.lo) == 0 || cmp(mid, window.hi) == 0)
    throw WindowIndivisible(
        "value-window holds a single key value with more pairs than the cap: "
        "raise the per-window pair cap");
  detail::window_descend(sq, KeyWindow{window.lo, mid}, max_pairs, scratch, f);
  detail::window_descend(sq, KeyWindow{mid, window.hi}, max_pairs, scratch, f);
}

} // namespace detail

// Calls f(const KeyWindow&, const WindowRanges&) for every non-empty window of
// a partition of the whole key range into windows of at most `max_pairs` pairs.
//
// `sq` must be sorted ascending -- checked, because an unsorted input does not
// fail, it silently drops pairs, and the caller is usually a sort whose
// comparator could be wrong. Windows arrive in ascending key order, which makes
// the whole pass deterministic: same input, same window sequence, same hits, in
// the same order.
//
// `scratch` is reused across every window and across calls; nothing here
// allocates.
template <typename F>
void for_each_window(std::span<const U256> sq, std::uint64_t max_pairs,
                     WindowRanges &scratch, F &&f) {
  const std::size_t d = sq.size();
  if (d < 2)
    return; // fewer than two points cannot make a pair
  if (d > 0xffffffffu)
    throw std::invalid_argument("for_each_window: D exceeds a uint32 index");
  if (max_pairs == 0)
    throw std::invalid_argument("for_each_window: max_pairs must be positive");
  for (std::size_t i = 1; i < d; ++i) {
    if (cmp(sq[i - 1], sq[i]) > 0)
      throw std::invalid_argument(
          "for_each_window: sq is not sorted ascending");
  }

  // The root spans every attainable key: the smallest pair sum is the two
  // smallest squares, the largest is the two largest, and the +1 makes the
  // half-open interval include it.
  const KeyWindow root{add(sq[0], sq[1]),
                       add(add(sq[d - 2], sq[d - 1]), u256_from(1))};
  detail::window_descend(sq, root, max_pairs, scratch, f);
}

} // namespace pte::core
