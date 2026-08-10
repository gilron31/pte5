// Host side of the champion meet: sort, plan, launch per window, confirm.
//
// The session is the object REVIEW.md A2(b)/A4(a) asked for and bench/gpu P7
// measured the need for: every device allocation -- sq, the hash table, the
// per-window ranges, the candidate buffer, the pinned staging block -- is made
// once and reused across every window. A champion run is thousands of launches,
// so a per-launch cudaMalloc or cudaHostAlloc would dominate it (P7: 9.2x on
// the staging path alone).
#include "gpu/champion/champion.h"

#include "core/window.h"
#include "gpu/champion/champion_impl.h"
#include "gpu/core/device_append_buffer.h"
#include "gpu/core/device_buffer.h"
#include "gpu/core/pinned_pool.h"
#include "gpu/core/types.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace pte::gpu {
namespace {

// Everything the run owns for its whole duration. Constructed once per
// champion_meet call and destroyed at the end of it; the buffers inside are
// what must not be per *window*.
class ChampionSession {
public:
  ChampionSession(std::size_t points, std::uint64_t slots,
                  std::size_t max_candidates)
      : sq_(points), lo_(points), hi_(points), table_(slots), overflow_(1),
        candidates_(static_cast<index_t>(max_candidates)) {}

  DeviceBuffer<core::U256> sq_;
  DeviceBuffer<std::uint32_t> lo_, hi_;
  DeviceBuffer<count_t> table_;
  DeviceBuffer<count_t> overflow_;
  DeviceAppendBuffer<ChampionCandidate> candidates_;
  PinnedPool pool_;
};

} // namespace

std::vector<core::PairCollision> champion_meet(std::span<const core::U128> nu,
                                               const ChampionParams &params,
                                               ChampionStats *stats) {
  std::vector<core::PairCollision> hits;
  const std::size_t d = nu.size();
  if (stats != nullptr)
    *stats = ChampionStats{static_cast<std::uint64_t>(d), 0, 0, 0, 0, 0};
  if (d < 4)
    return hits; // a collision is four distinct points

  if (d > kChampionMaxPoints)
    throw std::invalid_argument(
        "champion_meet: at most " + std::to_string(kChampionMaxPoints) +
        " points; a pair of indices shares one 64-bit table slot with the "
        "fingerprint");
  if (params.table_log < 4 || params.table_log > 40)
    throw std::invalid_argument("champion_meet: table_log out of range");
  if (!(params.max_load > 0.0 && params.max_load < 1.0))
    throw std::invalid_argument("champion_meet: max_load must be in (0,1)");

  const std::uint64_t slots = 1ull << params.table_log;
  const auto max_pairs =
      static_cast<std::uint64_t>(static_cast<double>(slots) * params.max_load);
  if (max_pairs == 0)
    throw std::invalid_argument("champion_meet: table too small for max_load");

  // Sorted by nu, which is sorted by nu^2 as well (nu >= 0 on canonized
  // points). `order` carries the caller's indexing through to the output.
  //
  // Written here rather than shared with cpu/champion on purpose: pte_gpu must
  // not link pte_cpu (REVIEW.md D5), and the two halves being independently
  // written is what makes their agreement evidence.
  std::vector<std::uint32_t> order(d);
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) {
    return core::cmp(nu[x], nu[y]) < 0;
  });
  for (std::size_t k = 1; k < d; ++k) {
    if (nu[order[k - 1]] == nu[order[k]])
      throw std::invalid_argument(
          "champion_meet: nu values must be distinct (a canonized point set "
          "has injective nu)");
  }

  std::vector<core::U256> sq(d);
  for (std::size_t k = 0; k < d; ++k)
    sq[k] = core::square(nu[order[k]]);

  ChampionSession session(d, slots, params.max_candidates);
  session.pool_.StageAndCopy(sq.data(), d, session.sq_);

  std::vector<ChampionCandidate> host_candidates(params.max_candidates);

  // The window planner writes its ranges straight into pinned memory -- the
  // producer filling the staging buffer rather than filling a vector and
  // copying it in (REVIEW.md A2(a)). Both ranges come from one Acquire because
  // the pool keeps one region live at a time.
  const std::span<std::uint32_t> staging = session.pool_.Acquire<std::uint32_t>(2 * d);
  core::WindowRanges ranges{staging.first(d), staging.subspan(d), 0};

  core::for_each_window(
      std::span<const core::U256>(sq), max_pairs, ranges,
      [&](const core::KeyWindow &, const core::WindowRanges &r) {
        if (stats != nullptr) {
          ++stats->windows;
          stats->pairs += r.pairs;
        }
        session.lo_.CopyFromHost(r.lo.data(), d);
        session.hi_.CopyFromHost(r.hi.data(), d);
        session.table_.Zero();
        session.overflow_.Zero();
        session.candidates_.Reset();

        LaunchChampionWindow(session.sq_, session.lo_, session.hi_,
                             session.table_, session.overflow_,
                             session.candidates_.view());

        // Reading the counter synchronizes, which is where an execution error
        // surfaces.
        const count_t found = session.candidates_.total();
        count_t overflow = 0;
        session.overflow_.CopyToHost(&overflow);
        if (overflow != 0)
          throw std::runtime_error(
              "champion_meet: hash table saturated (" +
              std::to_string(overflow) +
              " probe sequences exhausted); the window's result would be "
              "incomplete. Raise table_log or lower max_load");
        if (found > static_cast<count_t>(session.candidates_.capacity()))
          throw std::runtime_error(
              "champion_meet: " + std::to_string(found) +
              " candidates in one window exceeds max_candidates=" +
              std::to_string(params.max_candidates));

        if (found == 0)
          return;
        session.candidates_.CopyToHost(host_candidates.data(),
                                       static_cast<index_t>(found));
        if (stats != nullptr)
          stats->candidates += found;

        // Layer 2: the device compared a signature, this compares the key.
        for (count_t c = 0; c < found; ++c) {
          const ChampionCandidate &cand = host_candidates[static_cast<std::size_t>(c)];
          const core::U256 left = core::add(sq[cand.a], sq[cand.b]);
          const core::U256 right = core::add(sq[cand.c], sq[cand.d]);
          if (!(left == right)) {
            if (stats != nullptr)
              ++stats->false_positives;
            continue;
          }
          // Unreachable while nu is injective (checked above): equal keys over
          // distinct nu cannot share a point. Excluded anyway, because a
          // degenerate quadruple is not a solution.
          if (cand.a == cand.c || cand.a == cand.d || cand.b == cand.c ||
              cand.b == cand.d)
            continue;
          hits.push_back(core::PairCollision{order[cand.a], order[cand.b],
                                             order[cand.c], order[cand.d]}
                             .canonical());
        }
      });

  // Atomic append makes the device's hit order nondeterministic, so the output
  // is sorted before it is returned: a test that compared raw order against the
  // CPU reference would flake rather than fail.
  std::sort(hits.begin(), hits.end());
  if (stats != nullptr)
    stats->hits = hits.size();
  return hits;
}

} // namespace pte::gpu
