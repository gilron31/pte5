// The champion meet kernel: one value window, by open-addressing hash insert.
//
// ---------------------------------------------------------------------------
// The slot, and why it is one word
// ---------------------------------------------------------------------------
//
// A slot is a single 64-bit word:
//
//     [ occupied : 1 ][ fingerprint : 23 ][ i : 20 ][ j : 20 ]
//
// so the pair that owns a slot is stored *in* the slot. That is what makes the
// whole thing a single atomicCAS with no second step, and it is worth spelling
// out why the obvious alternative is worse.
//
// The obvious design stores the key in the slot and the pair beside it. Then a
// thread whose CAS reports "this key is already here" must read the partner's
// indices from an array the inserting thread has not necessarily written yet --
// a race with no atomic to close it, whose usual fix is a spin on a sentinel.
// Spinning inside a warp that also contains the writer is exactly the shape of
// a GPU deadlock, and it is the kind that appears at scale and never in a test.
//
// Packing the payload into the CAS word removes the second step entirely: if
// the CAS fails against an occupied slot, everything needed is already in the
// value it returned.
//
// **The consequence is that the comparison is on a signature, not on the key.**
// Two pairs are called a candidate when they share a slot sequence and 23 bits
// of fingerprint. That is sound in the direction that matters -- equal keys are
// identical 256-bit words, hence identical hashes, hence the same slot and the
// same fingerprint, so **no collision is ever missed** -- and unsound in the
// other, which the host's exact re-check absorbs.
//
// False-positive rate, since "23 bits" invites the wrong estimate: two keys are
// only ever compared if they land in the same probe sequence, so the effective
// width is table_log + 23 bits, not 23. At the sized window (2^26 slots, 3.3e7
// pairs) that is 49 bits against 2^49.6 pairs of pairs, i.e. **~1.6 false
// positives per window** -- rejected on the host at a cost of one 256-bit
// compare each.
//
// ---------------------------------------------------------------------------
// Concurrency
// ---------------------------------------------------------------------------
//
// Insertion is a linear probe from the home slot, and slots are never freed, so
// occupancy along a probe sequence is monotone. That gives the property the
// whole search rests on: **of two pairs with equal keys, exactly one sees the
// other.** The one that CASes later must have walked past the earlier one's
// slot to reach its own, and it only walks past a slot after observing it
// occupied -- so it observed the earlier entry and reported it. The one that
// CASed first never scans that far. Exactly one report, no duplicates, nothing
// missed.
//
// Spelled out, since "exactly one" is the part that makes the device answer
// equal to the reference answer rather than merely similar to it. Say A's CAS
// succeeded at slot p and B (equal key, hence the same probe sequence) later
// inserted at q. B advanced past p, and it advances past a slot only after the
// CAS there reported that slot occupied. Slots are never freed, so if B saw p
// occupied at time t then p was occupied for all later times -- and A's CAS at
// p succeeded, so A's write happened *before* t. B therefore compared against
// A's word and reported it. A, having stopped at p, never scanned as far as q
// and reports nothing. One report per unordered pair.
//
// One more detail makes the count match: a matching fingerprint reports and
// then keeps probing, because a third equal key has to find both of the first
// two. A group of g equal keys yields C(g,2) reports -- exactly what the CPU
// reference emits from an equal-key run.
#include "gpu/champion/champion_impl.h"

#include "core/uint.h"
#include "gpu/core/device_append_buffer.h"
#include "gpu/core/launch.cuh"

namespace pte::gpu {
namespace {

using core::U256;

constexpr int kIndexBits = kChampionIndexBits;
constexpr std::uint64_t kIndexMask = (1ull << kIndexBits) - 1;
constexpr int kFingerprintBits = 64 - 1 - 2 * kIndexBits;
constexpr std::uint64_t kOccupied = 1ull << 63;

static_assert(kFingerprintBits == 23, "the slot layout is 1 + 23 + 20 + 20");

// Probe limit. Not a correctness bound -- exhausting it is reported and treated
// as an error -- but a runtime bound, and it has to be far enough out that it
// never trips by chance. At load factor 0.5 the probability of a linear-probe
// cluster of length k decays like exp(-k(a - 1 - ln a)) = exp(-0.193 k), so
// clusters of 100 appear a few times in 3e7 inserts and clusters of 4096 have
// probability ~1e-343. It bounds the work; it does not bound the answer.
constexpr int kMaxProbes = 4096;

// Any deterministic function of the whole key would do; this is a 64-bit
// avalanche over all eight limbs. What matters is only that it depends on
// *every* limb -- a hash that ignored the high limbs would still be correct
// (equal keys still collide) but would put unrelated keys in the same probe
// sequence and turn the false-positive rate from 1.6 per window into millions.
__device__ inline std::uint64_t HashKey(const U256 &key) {
  std::uint64_t h = 0xcbf29ce484222325ull;
#pragma unroll
  for (int k = 0; k < 8; ++k) {
    h ^= key.w[k];
    h *= 0x100000001b3ull;
  }
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdull;
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53ull;
  h ^= h >> 33;
  return h;
}

template <int kThreads>
__global__ __launch_bounds__(kThreads) void ChampionWindowKernel(
    const U256 *__restrict__ sq, const std::uint32_t *__restrict__ lo,
    const std::uint32_t *__restrict__ hi, std::uint32_t num_points,
    count_t *__restrict__ table, std::uint64_t slot_mask,
    count_t *__restrict__ overflow, DeviceAppendView<ChampionCandidate> out) {
  const unsigned global_thread = blockIdx.x * kThreads + threadIdx.x;
  const std::uint32_t i = global_thread / 32u;
  const unsigned lane = global_thread % 32u;
  if (i >= num_points)
    return;

  const std::uint32_t begin = lo[i];
  const std::uint32_t end = hi[i];
  if (begin >= end)
    return; // most i are empty in most windows; two loads and out

  const U256 si = sq[i];
  const std::uint64_t packed_i = static_cast<std::uint64_t>(i) << kIndexBits;

  for (std::uint32_t j = begin + lane; j < end; j += 32u) {
    // Exact, because the window boundaries were decided on exact values: a
    // truncated key here would still hash consistently, but it would no longer
    // be the quantity the host confirms.
    const U256 key = core::add(si, sq[j]);
    const std::uint64_t h = HashKey(key);
    const std::uint64_t word =
        kOccupied | ((h >> (64 - kFingerprintBits)) << (2 * kIndexBits)) |
        packed_i | j;
    const std::uint64_t fingerprint_and_flag = word >> (2 * kIndexBits);

    std::uint64_t slot = h & slot_mask;
    int probe = 0;
    for (; probe < kMaxProbes; ++probe) {
      // The CAS *is* the probe. There is no load-then-insert gap to race on:
      // a failed CAS hands back the word that beat it, which is exactly what
      // the comparison below needs, and a slot once occupied stays occupied.
      const count_t prev =
          atomicCAS(table + slot, count_t{0}, static_cast<count_t>(word));
      if (prev == 0)
        break; // inserted; this pair is now the one a later equal key finds
      if ((prev >> (2 * kIndexBits)) == fingerprint_and_flag) {
        out.push(ChampionCandidate{
            i, j, static_cast<std::uint32_t>((prev >> kIndexBits) & kIndexMask),
            static_cast<std::uint32_t>(prev & kIndexMask)});
        // Report and keep going: a third equal key must find both of the first
        // two, and this pair still has to be inserted somewhere.
      }
      slot = (slot + 1) & slot_mask;
    }
    if (probe == kMaxProbes)
      atomicAdd(overflow, count_t{1});
  }
}

} // namespace

void LaunchChampionWindow(const DeviceBuffer<core::U256> &d_sq,
                          const DeviceBuffer<std::uint32_t> &d_lo,
                          const DeviceBuffer<std::uint32_t> &d_hi,
                          DeviceBuffer<count_t> &d_table,
                          DeviceBuffer<count_t> &d_overflow,
                          DeviceAppendView<ChampionCandidate> out) {
  const auto num_points = static_cast<std::uint32_t>(d_sq.count());
  if (num_points == 0 || d_table.count() == 0)
    return;
  const std::uint64_t slots = d_table.count();
  // One warp per point, so the grid is sized in warps rather than in threads.
  const unsigned long long warps = num_points;
  const unsigned blocks = static_cast<unsigned>(
      (warps * 32 + kChampionThreadsPerBlock - 1) / kChampionThreadsPerBlock);
  PTE_CUDA_LAUNCH_KERNEL((ChampionWindowKernel<kChampionThreadsPerBlock>),
                         blocks, kChampionThreadsPerBlock, d_sq.get(),
                         d_lo.get(), d_hi.get(), num_points, d_table.get(),
                         slots - 1, d_overflow.get(), out);
}

} // namespace pte::gpu
