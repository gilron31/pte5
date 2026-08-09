// Tests for the audit tooling.
//
// The audit is what stands in for an oracle at a bound too large to oracle, so
// it has to be checked itself: a verifier that always passes is worse than no
// verifier, because it looks like evidence.
#include "cpu/audit/audit.h"

#include "cpu/factorbase/factorbase.h"
#include "cpu/generate/generate.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

namespace {

using pte::cpu::Generator;
using pte::cpu::KnownCheck;
using pte::cpu::known_solutions;
using pte::cpu::norm_t;
using pte::cpu::SegmentCheck;

TEST(Audit, KnownSolutionsAreAscendingAndAttributed) {
  const auto known = known_solutions();
  ASSERT_GE(known.size(), 10u);
  for (std::size_t i = 1; i < known.size(); ++i)
    EXPECT_LT(known[i - 1].l1, known[i].l1) << "at index " << i;
  for (const auto &k : known) {
    EXPECT_NE(k.source, nullptr);
    EXPECT_GT(k.l1, 0);
  }
  // The two above 10^10 are the ones that make this table worth having: they
  // are unreachable by any bound this project has exhausted, so they are the
  // only evidence a 10^12 run would have about its own large-norm behaviour.
  const auto above = std::count_if(known.begin(), known.end(), [](const auto &k) {
    return k.l1 > 10'000'000'000LL;
  });
  EXPECT_EQ(above, 2) << "the exemplar and Bremner F2 k=7";
}

// Every known solution below 10^6 must be recovered. This is the audit working
// on a range where the answer is independently established.
TEST(Audit, RecoversEveryKnownSolutionInRange) {
  const KnownCheck c = pte::cpu::verify_known_solutions(1'000'000);
  EXPECT_EQ(c.in_range, 2u) << "67405 and 600445";
  EXPECT_EQ(c.found, c.in_range);
  EXPECT_TRUE(c.ok());
}

// And at a bound past every solution the exhaustive campaign confirmed.
TEST(Audit, RecoversEveryKnownSolutionToTenBillion) {
  const KnownCheck c = pte::cpu::verify_known_solutions(10'000'000'000LL);
  EXPECT_EQ(c.in_range, 8u) << "the eight the B = 10^10 campaign catalogued";
  EXPECT_EQ(c.found, 8u);
  EXPECT_TRUE(c.ok()) << "missing " << c.missing.size();
}

// The two above 10^10 are checkable without exhausting anything, which is the
// whole point: at B = 10^12 this is O(11) work, not O(2e10).
TEST(Audit, RecoversTheKnownSolutionsAboveTenBillion) {
  const KnownCheck c = pte::cpu::verify_known_solutions(100'000'000'000LL);
  EXPECT_EQ(c.in_range, 10u);
  EXPECT_EQ(c.found, 10u) << "33065996717 and 43937864165 must both certify";
  EXPECT_TRUE(c.ok());
}

// A bound below everything must report nothing rather than vacuously passing.
TEST(Audit, ReportsNothingInRangeBelowTheSmallestKnownSolution) {
  const KnownCheck c = pte::cpu::verify_known_solutions(1000);
  EXPECT_EQ(c.in_range, 0u);
  EXPECT_TRUE(c.ok());
}

TEST(Audit, SegmentReverificationAgreesWithTheReferencePath) {
  const Generator g(1'000'000);
  const SegmentCheck c = pte::cpu::verify_segments(g, 200);
  EXPECT_EQ(c.segments, 200u);
  EXPECT_GT(c.norms, 0u) << "a sample that walked no norms proves nothing";
  EXPECT_TRUE(c.ok()) << (c.disagreements.empty() ? "" : c.disagreements[0]);
}

// Deterministic sampling: an audit someone else cannot reproduce is not
// evidence. Same seed, same segments; different seed, different segments.
TEST(Audit, SamplingIsDeterministicAndSeedDependent) {
  const Generator g(1'000'000);
  const SegmentCheck a = pte::cpu::verify_segments(g, 64, 1);
  const SegmentCheck b = pte::cpu::verify_segments(g, 64, 1);
  const SegmentCheck c = pte::cpu::verify_segments(g, 64, 2);
  EXPECT_EQ(a.norms, b.norms);
  EXPECT_EQ(a.segments, b.segments);
  EXPECT_NE(a.norms, c.norms) << "two seeds picking identical work would make "
                                 "the seed parameter a lie";
}

} // namespace
