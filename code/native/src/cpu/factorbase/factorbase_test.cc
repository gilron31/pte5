// Tests for the segmented factor base.
//
// The load-bearing fixture is SegmentedAgreesWithMonolithic: the segmented
// sieve replaces the one every existing norm count was produced with, so
// "identical output" is what makes that replacement safe. Everything else in
// this file guards the artifact path, where a silently corrupt input would
// feed wrong primes into the enumeration.
#include "cpu/factorbase/factorbase.h"

#include "cpu/enumerate/enumerate.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using pte::cpu::factor_t;
using pte::cpu::FactorBase;
using pte::cpu::Gauss;
using pte::cpu::SegmentedSieve;

std::vector<factor_t> Sieved(factor_t p_max, std::size_t segment = 1u << 20) {
  std::vector<factor_t> out;
  SegmentedSieve(p_max, segment).for_each_split_prime([&](factor_t p) {
    out.push_back(p);
  });
  return out;
}

// A scratch path that cleans itself up even if a test fails.
class TempFile {
public:
  explicit TempFile(const char *name)
      : path_((std::filesystem::temp_directory_path() / name).string()) {}
  ~TempFile() { std::remove(path_.c_str()); }
  const std::string &path() const { return path_; }

private:
  std::string path_;
};

TEST(FactorBase, SegmentedAgreesWithMonolithic) {
  // The monolithic sieve is the one that produced every admissible-norm count
  // in the project. Replacing it is only safe if the output is identical.
  for (const factor_t limit : {factor_t{100}, factor_t{1000}, factor_t{65536},
                               factor_t{1000000}, factor_t{10000000}}) {
    EXPECT_EQ(Sieved(limit), pte::cpu::split_primes_up_to(limit))
        << "limit = " << limit;
  }
}

TEST(FactorBase, SegmentBoundariesDoNotChangeTheResult) {
  // Primes straddling a window edge are where an off-by-one lives. Widths that
  // are not powers of two, and one smaller than sqrt(p_max), exercise it.
  const auto reference = Sieved(200000);
  for (const std::size_t segment :
       {std::size_t{64}, std::size_t{97}, std::size_t{1024}, std::size_t{4096},
        std::size_t{100000}}) {
    EXPECT_EQ(Sieved(200000, segment), reference) << "segment = " << segment;
  }
}

TEST(FactorBase, HandlesDegenerateBounds) {
  EXPECT_TRUE(Sieved(0).empty());
  EXPECT_TRUE(Sieved(4).empty());
  EXPECT_EQ(Sieved(5), std::vector<factor_t>{5});
  EXPECT_THROW(SegmentedSieve(1000, 0), std::invalid_argument);
}

TEST(FactorBase, BuildDecomposesEveryPrime) {
  const FactorBase fb = FactorBase::Build(100000);
  ASSERT_FALSE(fb.empty());
  EXPECT_EQ(fb.primes(), pte::cpu::split_primes_up_to(100000));
  for (std::size_t i = 0; i < fb.size(); ++i) {
    const Gauss &z = fb.point(static_cast<pte::cpu::factor_index_t>(i));
    EXPECT_EQ(z.norm(), static_cast<std::int64_t>(fb.prime(
                            static_cast<pte::cpu::factor_index_t>(i))));
    EXPECT_GE(z.re, z.im);
    EXPECT_GT(z.im, 0);
  }
}

TEST(FactorBase, ArtifactRoundTrips) {
  const TempFile tmp("pte_fb_roundtrip.bin");
  const std::size_t written = FactorBase::WriteArtifact(tmp.path(), 1000000);

  const FactorBase built = FactorBase::Build(1000000);
  EXPECT_EQ(written, built.size());

  const FactorBase loaded = FactorBase::Load(tmp.path());
  ASSERT_EQ(loaded.size(), built.size());
  EXPECT_EQ(loaded.primes(), built.primes());
  for (std::size_t i = 0; i < loaded.size(); ++i) {
    const auto j = static_cast<pte::cpu::factor_index_t>(i);
    EXPECT_EQ(loaded.point(j), built.point(j)) << "at index " << i;
  }
}

TEST(FactorBase, LoadTruncatesToTheRequestedBound) {
  const TempFile tmp("pte_fb_truncate.bin");
  FactorBase::WriteArtifact(tmp.path(), 1000000);

  const FactorBase loaded = FactorBase::Load(tmp.path(), 100000);
  EXPECT_EQ(loaded.primes(), pte::cpu::split_primes_up_to(100000));

  // Asking for more than the file holds must fail loudly rather than quietly
  // enumerating a short base -- that would silently lose every norm whose
  // largest prime factor lies above the file's reach.
  EXPECT_THROW(FactorBase::Load(tmp.path(), 100000000), std::runtime_error);
}

TEST(FactorBase, LoadRejectsCorruptArtifacts) {
  // Wrong size: not a whole number of records.
  {
    const TempFile tmp("pte_fb_ragged.bin");
    std::FILE *f = std::fopen(tmp.path().c_str(), "wb");
    ASSERT_NE(f, nullptr);
    const char junk[12] = {};
    std::fwrite(junk, 1, sizeof(junk), f);
    std::fclose(f);
    EXPECT_THROW(FactorBase::Load(tmp.path()), std::runtime_error);
  }
  // Descending order, i.e. records swapped.
  {
    const TempFile tmp("pte_fb_unsorted.bin");
    std::FILE *f = std::fopen(tmp.path().c_str(), "wb");
    ASSERT_NE(f, nullptr);
    const std::uint32_t records[] = {3, 2, 2, 1}; // 13 then 5
    std::fwrite(records, sizeof(std::uint32_t), 4, f);
    std::fclose(f);
    EXPECT_THROW(FactorBase::Load(tmp.path()), std::runtime_error);
  }
  // a < b violates the canonical form.
  {
    const TempFile tmp("pte_fb_uncanonical.bin");
    std::FILE *f = std::fopen(tmp.path().c_str(), "wb");
    ASSERT_NE(f, nullptr);
    const std::uint32_t records[] = {1, 2};
    std::fwrite(records, sizeof(std::uint32_t), 2, f);
    std::fclose(f);
    EXPECT_THROW(FactorBase::Load(tmp.path()), std::runtime_error);
  }
}

// The shipped artifact is an input to the campaign, so its first records are
// worth pinning: a regression in WriteArtifact that changed the layout would
// otherwise only surface as wrong answers much later.
TEST(FactorBase, ShippedArtifactPrefixMatchesAFreshBuild) {
  const std::filesystem::path shipped =
      std::filesystem::path(PTE_DATA_DIR) / "factorbase_1e12.bin";
  if (!std::filesystem::exists(shipped))
    GTEST_SKIP() << "data/factorbase_1e12.bin not present";

  const FactorBase prefix = FactorBase::Load(shipped.string(), 1000000);
  const FactorBase fresh = FactorBase::Build(1000000);
  ASSERT_EQ(prefix.size(), fresh.size());
  EXPECT_EQ(prefix.primes(), fresh.primes());
  for (std::size_t i = 0; i < prefix.size(); ++i) {
    const auto j = static_cast<pte::cpu::factor_index_t>(i);
    EXPECT_EQ(prefix.point(j), fresh.point(j)) << "at index " << i;
  }
}

} // namespace
