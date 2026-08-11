// G1, the local model: depth-d towers over Z/p^k, with a cap on node valuations.
//
// ---------------------------------------------------------------------------
// Why mod p^k, and what the cap is for
// ---------------------------------------------------------------------------
//
// drafts/gpu_implementation_plan.md 2.2 rejects two candidate local models and
// leaves sigma_p undefined:
//
//   distinct roots     counts towers whose 2^d roots stay distinct mod p. Gives
//                      zero depth-5 towers at 25 of the 29 primes <= 113, which
//                      read naively is an unconditional nonexistence proof. It
//                      is not: the 16 roots of the section-7 exemplar E_4 are
//                      distinct over Z and collide mod p at 16 of those primes.
//   with multiplicity  counts towers that split mod p at all. Once the node set
//                      collapses to {0} the only completion is L_d = 0, i.e.
//                      E_d = E_{d-1}^2 -- a tower, but not a solution.
//
// **They are not two competing answers.** They are the two extreme cases of one
// computation, and this file is that computation.
//
// A solution is a Z-point, hence a Z_p-point, of the *nondegenerate* locus, and
// nondegenerate means every node value is nonzero in Z_p -- that its valuation
// is **finite**. Mod p you can see "zero" and "nonzero"; you cannot see
// "finite". So bound the valuation instead: fix a cap v and admit only node
// values t with v_p(t) <= v, under the convention v_p(0) = k, which is the
// largest valuation Z/p^k can tell apart from infinity.
//
// At k = 1 the two rejected models are exactly the two extreme caps:
//
//   max_valuation = 0    v_p(t) = 0 means t != 0 mod p    the distinct-root model
//   max_valuation = k    no restriction at all            the multiplicity model
//
// and tower_test.cc asserts precisely that, against the numbers
// drafts/verify_gpu_plan.py prints. **The cap is the missing parameter and k is
// the missing precision**; drafts/local_model_plan.md is the argument.
//
// ---------------------------------------------------------------------------
// What it computes
// ---------------------------------------------------------------------------
//
// The recursion of paper.md 2.1, unchanged except that the ring is Z/p^k rather
// than F_p:
//
//     T_0 = Z/p^k,   A_j = { t in T_j : -t in T_j, v_p(t) <= cap },
//     T_{j+1} = { t^2 - L_{j+1} : t in A_j },
//
// with the depth-d tower splitting iff 0 in T_d. The elements of A_j are the
// node values at tree level j+1, which is why capping A is capping every node.
//
// Both branches of drafts/local_model_plan.md 1 read off this one run:
//
//   obstruction   V_p(k,d) = min over caps admitting a tower. Growing without
//                 bound in k means every Z_p-tower is degenerate, hence no E_d
//                 over Z. Stabilising at v_0 means every E_d has a node value
//                 divisible by p^{v_0} -- paper.md 4's ladder, computed.
//   density       towers[d-1] is N_k, the point count mod p^k that feeds the
//                 Poincare series.
//
// ---------------------------------------------------------------------------
// Cost, and the one number that decides how far this goes
// ---------------------------------------------------------------------------
//
// A state is a subset T of Z/m as an m-bit mask, m = p^k. The transition
// T' = { s - L : s in S } is a **cyclic rotation of the S-mask by L**, so one
// level costs states * m * m/64 word-operations and the memory is
// states * m/8 bytes.
//
// Everything therefore hinges on how `states` grows. Measured at p = 3, it is
// O(p^k) and not O(p^{2k}) -- and the reason is structural: every state at
// level >= 1 has the form S - L, so each level's state set is a **union of
// complete translation orbits**. Set `count_orbits` to measure it: `orbits` is
// the number of orbits and `orbit_misses` must be 0, which is the assertion
// that the structure is real. If it holds, the orbit representative -- not the
// state -- is the object a larger run should store, and drafts/'s 2.3 GB at
// 5^7 becomes kilobytes.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace pte::cpu {

// "No restriction on node valuations." Numerically any cap >= k does this,
// since v_p(0) is capped at k; this spelling says so at the call site.
inline constexpr std::uint32_t kNoValuationCap = 0xffffffffu;

struct TowerParams {
  std::uint32_t p = 3;
  std::uint32_t k = 1;
  std::uint32_t depth = 5;

  // Node values t are admitted only when v_p(t) <= max_valuation. Clamped to k.
  std::uint32_t max_valuation = kNoValuationCap;

  // Exact tuple counts. Off is not an optimization -- it is what lets a run
  // exceed the 64-bit accumulator, since the count is bounded by m^depth and
  // that is what `tower_counts` refuses to overflow silently. The obstruction
  // branch needs `splitting_states` only, so it can and should turn this off.
  bool count_towers = true;

  // The translation-orbit diagnostic described above. Costs an extra
  // states * m * m/64 pass, so it is opt-in.
  bool count_orbits = false;

  // Fail fast rather than exhaust memory: each state costs m/8 bytes.
  std::uint64_t max_states = 1u << 20;
};

struct TowerCounts {
  // towers[j] = the number of prefixes (L_1..L_{j+1}) in (Z/p^k)^{j+1} whose
  // depth-(j+1) tower splits. Empty when count_towers is false.
  std::vector<std::uint64_t> towers;

  // states[j] = distinct T-sets entering level j+1. states[0] is always 1
  // (the full set), so this is directly comparable with the "states/level"
  // column drafts/verify_gpu_plan.py prints.
  std::vector<std::uint64_t> states;

  // splitting_states[j] = how many of those states admit at least one L
  // completing a depth-(j+1) tower. Always computed, never overflows, and
  // > 0 is the whole of the existence question.
  std::vector<std::uint64_t> splitting_states;

  // Only when count_orbits: states[j] grouped under translation, and the count
  // of (state, L) pairs whose rotation is *not* itself a state. The latter is
  // expected to be 0 at every level; a nonzero entry falsifies the orbit
  // structure the header rests on.
  std::vector<std::uint64_t> orbits;
  std::vector<std::uint64_t> orbit_misses;

  bool splits() const {
    return !splitting_states.empty() && splitting_states.back() > 0;
  }
};

// Thrown when a run would exceed max_states, or when exact counting would
// overflow 64 bits. Both are sizing errors, not arithmetic that quietly wraps.
struct TowerTooLarge : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Throws std::invalid_argument on a non-prime p, k = 0, depth = 0, or a
// modulus past the sanity bound; TowerTooLarge on the two sizing failures.
TowerCounts tower_counts(const TowerParams &params);

struct ValuationFloor {
  std::uint32_t p = 0, k = 0, depth = 0;

  // V_p(k,d): the least cap admitting a depth-d tower over Z/p^k. Equal to k
  // means "nothing nondegenerate is visible at this precision" -- the value
  // that grows into a nonexistence proof if it keeps growing with k, and the
  // value that means nothing at all if k is below the precision a real solution
  // needs (drafts/local_model_plan.md 2 fixes that floor from paper.md 4).
  std::uint32_t floor = 0;

  // towers_by_cap[v] for v = 0..floor, so the run that finds the floor also
  // hands back the valuation distribution the density branch wants.
  std::vector<std::uint64_t> towers_by_cap;

  // From the run at cap = floor.
  std::vector<std::uint64_t> states;

  bool nondegenerate() const { return floor < k; }
};

// Ascending scan over caps with an early exit, which is correct because the
// admissible-value sets are nested and the whole recursion is monotone in them:
// a larger cap can only add towers.
//
// With base.count_towers off it first probes cap = k-1, which settles "nothing
// nondegenerate at this precision" in one run rather than k. That is the common
// answer at large k and it is the only reason the sweep reaches 3^8 at all.
ValuationFloor valuation_floor(std::uint32_t p, std::uint32_t k,
                               std::uint32_t depth,
                               const TowerParams &base = {});

} // namespace pte::cpu
