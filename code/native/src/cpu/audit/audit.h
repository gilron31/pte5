// Auditing an exhaustive run at a bound nobody can cross-check.
//
// Every strong correctness fixture in this project lives at B = 10^6 or 10^7,
// where an independently written oracle can be run over the whole space. At
// B = 10^12 there is no such oracle -- 2e10 norms, and the reference enumerator
// is 2.5x slower than the thing it would be checking. So the assurance has to
// change shape rather than scale.
//
// Three mechanisms here, in increasing order of what they cost:
//
//   1. KNOWN SOLUTIONS. An E_4 the run *must* find. These come from outside
//      this codebase entirely -- Bremner's parameterized families and the
//      sporadic solutions recorded in code/old -- so they are evidence about
//      the search, not about its self-consistency. Crucially, a known solution
//      can be checked at *any* bound by enumerating one norm rather than all of
//      them: verify_known_solutions is O(number of known solutions), not O(B).
//
//   2. SEGMENT RE-VERIFICATION. Re-run a sample of segments through the
//      *reference* path -- Enumerator's rebuild-from-scratch point sets and the
//      naive exact meet_oracle -- and require identical results. Segments are
//      independently enumerable, so this audits a real slice of a run that is
//      far too large to audit whole. This is what makes an exclusion claim
//      checkable by someone else.
//
//   3. NORM COUNT. The count of admissible norms is the one aggregate that a
//      silent miss must move. Comparing it against an independent count is the
//      cheapest evidence that the DFS is not dropping a class of norms.
//
// None of these proves exhaustiveness. Together they make the failure modes
// that an exclusion claim dies of -- a dropped prime, a lost non-squarefree
// norm, a truncated factor base -- observable rather than silent.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/generate/generate.h"

#include <cstddef>
#include <string>
#include <vector>

namespace pte::cpu {

// An E_4 known from outside this codebase.
struct KnownSolution {
  norm_t l1;
  const char *source;
};

// Every E_4 this project knows of from an independent source, ascending by L1.
//
// Provenance matters more than completeness here, so it is recorded per entry:
//
//   - "bremner F1 k=N" / "bremner F2 k=N" come from Bremner's two parameterized
//     families, evaluated in code/old (utils.bremner_F1_point / _F2_point).
//     Being families, they generate solutions at *arbitrary* size -- which is
//     what makes them usable above any bound this project can exhaust.
//   - "bremner sporadic" are the five isolated solutions listed alongside those
//     families.
//   - "gilro sporadic" are from the earlier champion-norm enumeration.
//
// The two entries above 10^10 are the ones that matter for a design-bound run:
// 33,065,996,717 (the paper's exemplar) and 43,937,864,165 must both appear in
// an exhaustive catalog to 10^12, and their absence would be a bug that no
// self-consistency check could see.
std::vector<KnownSolution> known_solutions();

// Result of checking the known solutions against a bound.
struct KnownCheck {
  std::size_t in_range = 0;   // known solutions with L1 <= bound
  std::size_t found = 0;      // ...that the pipeline recovered and certified
  std::vector<norm_t> missing; // ...that it did not
  bool ok() const { return missing.empty(); }
};

// For each known solution at or below `bound`, enumerate *only* that norm,
// run the production meet, and require a certified primitive hit.
//
// Costs one point-set construction per known solution, so it is instant at any
// bound. It does not check that the campaign's DFS reaches that norm -- see
// verify_segments for that -- it checks that the pipeline recognizes the
// solution when handed the norm.
KnownCheck verify_known_solutions(norm_t bound);

// Result of re-running segments through the reference implementation.
struct SegmentCheck {
  std::size_t segments = 0;
  std::size_t norms = 0;
  std::size_t solutions = 0;
  std::vector<std::string> disagreements;
  bool ok() const { return disagreements.empty(); }
};

// Re-verify `count` segments, sampled deterministically from `seed`, against
// the reference path: point sets rebuilt by point_set_from, meet by the naive
// exact meet_oracle. Any difference in the solution set is a disagreement.
//
// Deterministic sampling so an audit is reproducible and citable -- "segments
// sampled with seed 1 at B = 10^12" is a statement someone else can re-run.
SegmentCheck verify_segments(const Generator &g, std::size_t count,
                             std::uint64_t seed = 1);

} // namespace pte::cpu
