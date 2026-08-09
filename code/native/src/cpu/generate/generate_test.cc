// Tests for the DFS spine.
//
// The load-bearing fixture is SpineMatchesRebuildOnEveryNorm: the spine builds
// point sets incrementally, the reference builds each from {1,0}, and they must
// produce identical vectors -- not merely sets of the same size. Everything the
// meet does rests on those points being canonized, distinct, and complete, and
// the incremental construction is where all three could quietly stop being true.
//
// Note what is *not* shared: Enumerator keeps its own factor base, its own
// largest-prime cutoff and its own point_set_from. It is the independently
// written reference, in the same way meet_oracle is for the meet.
#include "cpu/generate/generate.h"

#include "cpu/enumerate/enumerate.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/gaussian/gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using pte::cpu::factor_t;
using pte::cpu::FactorBase;
using pte::cpu::Gauss;
using pte::cpu::Generator;
using pte::cpu::norm_t;
using pte::cpu::NormView;
using pte::cpu::Spine;

// Every admissible norm the reference produces, with its point set.
std::map<norm_t, std::vector<Gauss>> ReferenceNorms(norm_t bound) {
  const pte::cpu::Enumerator e(bound);
  std::map<norm_t, std::vector<Gauss>> out;
  e.for_each_norm([&](const pte::cpu::DecomposedNorm &norm) {
    out.emplace(norm.l1, e.point_set(norm));
  });
  return out;
}

TEST(Spine, MatchesRebuildOnEveryNormToOneMillion) {
  const auto reference = ReferenceNorms(1'000'000);
  ASSERT_EQ(reference.size(), 11'083u)
      << "the admissible-norm count must not move";

  const Generator g(1'000'000);
  std::size_t seen = 0;
  g.for_each_norm([&](const NormView &norm) {
    const auto it = reference.find(norm.l1);
    ASSERT_NE(it, reference.end()) << "L1 = " << norm.l1 << " is not admissible";
    // Identical vectors, not just identical sizes: the spine sorts and dedupes
    // at every level, so element order is part of the contract the meet's
    // injectivity precondition is checked against.
    const std::vector<Gauss> got(norm.points.begin(), norm.points.end());
    EXPECT_EQ(got, it->second) << "L1 = " << norm.l1;

    // nu must be in point order and injective. The meet reads this span
    // directly and drops its disjointness test on the strength of injectivity,
    // so a permuted or collapsed nu would not fail loudly -- it would report
    // solutions whose four indices are not four distinct points.
    const std::vector<norm_t> want = pte::cpu::nu_values(it->second);
    const std::vector<norm_t> have(norm.nu.begin(), norm.nu.end());
    ASSERT_EQ(have, want) << "L1 = " << norm.l1;
    std::vector<norm_t> sorted = have;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_EQ(std::adjacent_find(sorted.begin(), sorted.end()), sorted.end())
        << "nu is not injective at L1 = " << norm.l1;
    ++seen;
  });
  EXPECT_EQ(seen, reference.size()) << "the spine must emit every norm once";
}

// Non-squarefree norms are where the incremental construction is hardest: the
// point set stops doubling and the deduplication becomes many-to-one. They are
// also what Prop. 5.1's exhaustiveness rests on, so losing them would be a
// silent hole in the paper rather than a slow program.
TEST(Spine, NonSquarefreeNormsSurviveTheIncrementalBuild) {
  const auto reference = ReferenceNorms(200'000);
  const Generator g(200'000);

  std::size_t repeated = 0;
  g.for_each_norm([&](const NormView &norm) {
    const std::vector<factor_t> primes = g.primes_of(norm);
    const bool has_repeat =
        std::adjacent_find(primes.begin(), primes.end()) != primes.end();
    if (!has_repeat)
      return;
    ++repeated;
    const std::vector<Gauss> got(norm.points.begin(), norm.points.end());
    EXPECT_EQ(got, reference.at(norm.l1)) << "L1 = " << norm.l1;
  });
  EXPECT_GT(repeated, 400u) << "about a fifth of admissible norms repeat a prime";
}

// The three-way collapse from the header, in isolation. U(25) has two points;
// multiplying by 2+i and its conjugate gives four products of which three
// canonize to the same point. A dedup that assumed "at most one duplicate"
// would return three points for a norm that has two.
TEST(Spine, DeduplicationIsManyToOneNotPairwise) {
  const FactorBase base = FactorBase::Build(100);
  Spine spine(base);

  const auto index_of = [&](factor_t p) {
    const auto &primes = base.primes_vector();
    return static_cast<pte::cpu::factor_index_t>(
        std::lower_bound(primes.begin(), primes.end(), p) - primes.begin());
  };
  const auto five = index_of(5);

  spine.push(five); // 5
  EXPECT_EQ(spine.value(), 5);
  EXPECT_EQ(spine.points().size(), 1u);

  spine.push(five); // 25: d = 3, |U| = 2
  EXPECT_EQ(spine.value(), 25);
  ASSERT_EQ(spine.points().size(), 2u);
  EXPECT_EQ(spine.points()[0], (Gauss{4, 3}));
  EXPECT_EQ(spine.points()[1], (Gauss{5, 0}));

  spine.push(five); // 125: four products, three of them (10,5)
  EXPECT_EQ(spine.value(), 125);
  ASSERT_EQ(spine.points().size(), 2u)
      << "the collapse is three-to-one, not a single duplicate pair";
  EXPECT_EQ(spine.points()[0], (Gauss{10, 5}));
  EXPECT_EQ(spine.points()[1], (Gauss{11, 2}));

  spine.pop();
  EXPECT_EQ(spine.points().size(), 2u) << "pop must restore the parent's set";
  EXPECT_EQ(spine.points()[1], (Gauss{5, 0}));
}

TEST(Spine, PopRestoresEveryAncestor) {
  const FactorBase base = FactorBase::Build(1000);
  Spine spine(base);

  // Record each level's points on the way down, then check them on the way up.
  std::vector<std::vector<Gauss>> going_down;
  for (pte::cpu::factor_index_t i = 0; i < 6; ++i) {
    spine.push(i);
    going_down.emplace_back(spine.points().begin(), spine.points().end());
  }
  for (int i = 5; i >= 0; --i) {
    const std::vector<Gauss> got(spine.points().begin(), spine.points().end());
    EXPECT_EQ(got, going_down[static_cast<std::size_t>(i)]) << "at depth " << i + 1;
    spine.pop();
  }
  EXPECT_EQ(spine.depth(), 0);
  EXPECT_EQ(spine.points().size(), 1u) << "U(1) = {1}";
  EXPECT_EQ(spine.points()[0], (Gauss{1, 0}));
}

// Depth is driven by the smallest split prime, so 5^k is the worst case: at
// B = 10^12 that is seventeen levels. A budget derived from chains of distinct
// primes -- an earlier draft used 13 -- is not merely tight, it is wrong.
TEST(Spine, DepthBudgetCoversTheSmallestPrimeChain) {
  const FactorBase base = FactorBase::Build(100);
  Spine spine(base);
  for (int k = 0; k < pte::cpu::kMaxSpineDepth; ++k)
    spine.push(0); // 5^k
  EXPECT_EQ(spine.depth(), pte::cpu::kMaxSpineDepth);

  // And the budget is enforced, not merely documented.
  EXPECT_THROW(spine.push(0), std::length_error);

  // 5^17 <= 10^12 < 5^18, so the design bound needs exactly 17.
  norm_t v = 1;
  int k = 0;
  while (v <= pte::cpu::kMaxNorm / 5) {
    v *= 5;
    ++k;
  }
  EXPECT_EQ(k, 17);
  EXPECT_LE(k, pte::cpu::kMaxSpineDepth);
}

TEST(Generator, HandlesBoundsBelowTheSmallestAdmissibleNorm) {
  std::size_t count = 0;
  Generator(1000).for_each_norm([&](const NormView &) { ++count; });
  EXPECT_EQ(count, 0u) << "5 * 13 * 17 = 1105 is the smallest admissible norm";
  EXPECT_THROW(Generator(-1), std::invalid_argument);
}

// A base that does not reach bound / 65 must be refused. This is the failure an
// exclusion claim cannot survive: the run completes, reports fewer norms, and
// looks exactly like a correct one.
TEST(Generator, RejectsAFactorBaseThatIsTooShort) {
  EXPECT_THROW(Generator(1'000'000, FactorBase::Build(1000)),
               std::invalid_argument);

  // Reaching exactly bound / 65 is enough...
  const Generator ok(1'000'000, FactorBase::Build(1'000'000 / 65));
  std::size_t count = 0;
  ok.for_each_norm([&](const NormView &) { ++count; });
  EXPECT_EQ(count, 11'083u);

  // ...and so is a base that overshoots, which must not change the answer or
  // cost the top-level loop a single extra prime.
  const Generator over(1'000'000, FactorBase::Build(1'000'000));
  std::size_t over_count = 0;
  over.for_each_norm([&](const NormView &) { ++over_count; });
  EXPECT_EQ(over_count, 11'083u);
}

TEST(Generator, DMatchesThePointCountFormula) {
  const Generator g(1'000'000);
  std::map<int, std::size_t> hist;
  g.for_each_norm([&](const NormView &norm) {
    // The bucket the meet dispatches on must agree with prod(b_i + 1) / 2,
    // computed independently from the factorization.
    ASSERT_EQ(static_cast<std::size_t>(norm.D()),
              pte::cpu::point_count(norm.factors))
        << "L1 = " << norm.l1;
    ASSERT_GE(norm.D(), 4) << "L1 = " << norm.l1;
    ++hist[norm.D()];
  });
  // The measured distribution, including the odd sizes the design doc's
  // "D <= 8" column dropped.
  EXPECT_EQ(hist[4], 8805u);
  EXPECT_EQ(hist[5], 154u);
  EXPECT_EQ(hist[6], 1344u);
  EXPECT_EQ(hist[7], 7u);
  EXPECT_EQ(hist[8], 644u);
  EXPECT_EQ(hist[16], 1u);
}

} // namespace
