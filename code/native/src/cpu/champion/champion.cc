#include "cpu/champion/champion.h"

#include "core/uint.h"
#include "core/window.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

namespace pte::cpu {
namespace {

// One pair, with its exact key. 40 bytes; this is what bounds the reference's
// memory to max_pairs_per_window entries.
struct Entry {
  core::U256 key;
  std::uint32_t i, j; // indices in *sorted* order, not the caller's
};

} // namespace

std::vector<core::PairCollision> champion_meet(std::span<const core::U128> nu,
                                               const ChampionParams &params,
                                               ChampionStats *stats) {
  std::vector<core::PairCollision> hits;
  const std::size_t d = nu.size();
  if (stats != nullptr)
    *stats = ChampionStats{static_cast<std::uint64_t>(d), 0, 0, 0};
  if (d < 4)
    return hits; // a collision is four distinct points

  // Sort by nu, which is the same order as by nu^2 (nu >= 0 on canonized
  // points), and remember where each point came from so hits can be reported in
  // the caller's indexing.
  std::vector<std::uint32_t> order(d);
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) {
    return core::cmp(nu[x], nu[y]) < 0;
  });

  // The injectivity precondition, checked once in O(D) now that the values are
  // adjacent. Cheap here, and the alternative is emitting quadruples that share
  // a point -- see the header.
  for (std::size_t k = 1; k < d; ++k) {
    if (nu[order[k - 1]] == nu[order[k]])
      throw std::invalid_argument(
          "champion_meet: nu values must be distinct (a canonized point set "
          "has injective nu); a repeated value means the points were not "
          "canonized or not deduplicated");
  }

  std::vector<core::U256> sq(d);
  for (std::size_t k = 0; k < d; ++k)
    sq[k] = core::square(nu[order[k]]);

  // Reserve for the smaller of the cap and the whole problem, so a four-point
  // test does not allocate the 40 MB the champion-scale default implies.
  std::vector<Entry> entries;
  entries.reserve(static_cast<std::size_t>(
      std::min<std::uint64_t>(params.max_pairs_per_window,
                              static_cast<std::uint64_t>(d) * (d - 1) / 2)));

  core::WindowScratch scratch(d);
  core::WindowRanges ranges = scratch.ranges();
  core::for_each_window(
      std::span<const core::U256>(sq), params.max_pairs_per_window, ranges,
      [&](const core::KeyWindow &, const core::WindowRanges &r) {
        if (stats != nullptr) {
          ++stats->windows;
          stats->pairs += r.pairs;
        }
        entries.clear();
        for (std::size_t i = 0; i < d; ++i) {
          const core::U256 &si = sq[i];
          for (std::uint32_t j = r.lo[i]; j < r.hi[i]; ++j)
            entries.push_back(
                Entry{core::add(si, sq[j]), static_cast<std::uint32_t>(i), j});
        }

        // Sort and read off equal-key runs. This is the same shape as
        // meet_hashed one level up, and for the same reason: it turns an
        // O(pairs^2) comparison into one pass over sorted keys. The tie-break
        // on (i, j) is not cosmetic -- it makes the output order a function of
        // the input alone, so two runs of the same champion agree exactly.
        std::sort(entries.begin(), entries.end(),
                  [](const Entry &x, const Entry &y) {
                    const int c = core::cmp(x.key, y.key);
                    if (c != 0)
                      return c < 0;
                    return x.i != y.i ? x.i < y.i : x.j < y.j;
                  });

        for (std::size_t a = 0; a < entries.size();) {
          std::size_t b = a + 1;
          while (b < entries.size() && entries[b].key == entries[a].key)
            ++b;
          // Runs of length one are the overwhelming majority and cost one
          // comparison; only a genuine collision reaches the inner loops.
          for (std::size_t p = a; p + 1 < b; ++p)
            for (std::size_t q = p + 1; q < b; ++q) {
              const Entry &u = entries[p];
              const Entry &v = entries[q];
              // Unreachable while nu is injective, which is checked above --
              // but a shared index would be a degenerate "solution" rather
              // than a wrong one, so it costs nothing to exclude it here too.
              if (u.i == v.i || u.i == v.j || u.j == v.i || u.j == v.j)
                continue;
              hits.push_back(core::PairCollision{order[u.i], order[u.j],
                                                 order[v.i], order[v.j]}
                                 .canonical());
            }
          a = b;
        }
      });

  // Sorted, but deliberately *not* deduplicated: every pair lands in exactly
  // one window (core/window.h) and every equal-key run emits each pair of pairs
  // once, so a duplicate here would be a bug in one of those two -- and a
  // std::unique that quietly swallowed it would hide the only symptom.
  std::sort(hits.begin(), hits.end());
  if (stats != nullptr)
    stats->hits = hits.size();
  return hits;
}

} // namespace pte::cpu
