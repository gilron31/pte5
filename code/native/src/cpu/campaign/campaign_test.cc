// Tests for the campaign.
//
// Two properties matter and neither is about speed:
//
//   1. The catalog does not depend on the thread count. If it did, "the
//      catalog" would not be a well-defined object and no two runs could be
//      compared.
//   2. A run that stops and resumes produces the same catalog as one that never
//      stopped. This is the property an exclusion claim rests on: a campaign
//      that lost a batch across a restart reports fewer solutions and looks
//      exactly like a correct one.
#include "cpu/campaign/campaign.h"

#include "cpu/confirm/confirm.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/generate/generate.h"
#include "cpu/meet/meet.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

using pte::cpu::Campaign;
using pte::cpu::CampaignOptions;
using pte::cpu::CampaignStats;

// A scratch directory that cleans itself up even when a test fails.
class TempDir {
public:
  TempDir() {
    path_ = std::filesystem::temp_directory_path() /
            ("pte_campaign_" + std::to_string(::getpid()) + "_" +
             std::to_string(counter_++));
    std::filesystem::create_directories(path_);
  }
  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
  }
  std::string file(const char *name) const { return (path_ / name).string(); }

private:
  std::filesystem::path path_;
  static inline int counter_ = 0;
};

std::vector<std::string> ReadLines(const std::string &path) {
  std::ifstream in(path);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line))
    lines.push_back(line);
  return lines;
}

CampaignOptions Base(const TempDir &dir, const char *out) {
  CampaignOptions o;
  o.bound = 1'000'000;
  o.out = dir.file(out);
  o.batch = 64;
  o.quiet = true;
  return o;
}

TEST(Campaign, FindsTheKnownSolutionsAndCountsEveryNorm) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "catalog.csv");
  o.threads = 4;
  const CampaignStats s = Campaign(o).run();

  EXPECT_TRUE(s.complete);
  EXPECT_EQ(s.norms, 11'083u) << "the admissible-norm count must not move";
  EXPECT_EQ(s.solutions, 2u) << "the known E_4 count below 10^6";
  EXPECT_EQ(s.uncertified, 0u) << "every catalog row must carry a proof";

  const std::vector<std::string> lines = ReadLines(o.out);
  ASSERT_EQ(lines.size(), 3u) << "header plus two solutions";
  EXPECT_EQ(lines[0], Campaign::header());
  // Certified flag is the last column.
  for (std::size_t i = 1; i < lines.size(); ++i)
    EXPECT_EQ(lines[i].substr(lines[i].rfind(',') + 1), "1") << lines[i];
  EXPECT_EQ(lines[1].substr(0, lines[1].find(',')), "67405");
  EXPECT_EQ(lines[2].substr(0, lines[2].find(',')), "600445");
}

// Segments are written back in segment order, so the byte stream must not
// depend on how many threads raced to produce it.
TEST(Campaign, CatalogIsIndependentOfThreadCount) {
  const TempDir dir;
  std::vector<std::string> reference;
  for (const int threads : {1, 2, 8}) {
    CampaignOptions o = Base(dir, "t.csv");
    o.out = dir.file(threads == 1   ? "t1.csv"
                     : threads == 2 ? "t2.csv"
                                    : "t8.csv");
    o.threads = threads;
    const CampaignStats s = Campaign(o).run();
    EXPECT_EQ(s.norms, 11'083u) << "threads = " << threads;
    EXPECT_EQ(s.solutions, 2u) << "threads = " << threads;

    const std::vector<std::string> lines = ReadLines(o.out);
    if (reference.empty())
      reference = lines;
    else
      EXPECT_EQ(lines, reference) << "threads = " << threads;
  }
}

// The load-bearing fixture. Stop after two batches, resume, and land on the
// same catalog as an uninterrupted run -- byte for byte, and with the same
// counters.
TEST(Campaign, ResumeIsEquivalentToAnUninterruptedRun) {
  const TempDir dir;

  CampaignOptions whole = Base(dir, "whole.csv");
  whole.threads = 4;
  const CampaignStats full = Campaign(whole).run();
  ASSERT_TRUE(full.complete);

  CampaignOptions part = Base(dir, "part.csv");
  part.threads = 4;
  part.stop_after_batches = 2;
  const CampaignStats first = Campaign(part).run();
  EXPECT_FALSE(first.complete);
  EXPECT_GT(first.segments_done, 0u);
  EXPECT_LT(first.segments_done, full.segments_done);

  // Resume, possibly several times, to prove the checkpoint is re-enterable
  // rather than usable once.
  CampaignStats last = first;
  int rounds = 0;
  while (!last.complete) {
    CampaignOptions more = part;
    more.resume = true;
    more.threads = 3; // a different width, to catch order dependence
    more.batch = 37;  // and different boundaries, which must not matter
    last = Campaign(more).run();
    ASSERT_LT(++rounds, 100) << "resume made no progress";
  }

  EXPECT_EQ(last.norms, full.norms);
  EXPECT_EQ(last.candidates, full.candidates);
  EXPECT_EQ(last.filtered, full.filtered);
  EXPECT_EQ(last.imprimitive, full.imprimitive);
  EXPECT_EQ(last.solutions, full.solutions);
  EXPECT_EQ(last.segments_done, full.segments_done);
  EXPECT_EQ(ReadLines(part.out), ReadLines(whole.out));
}

// Rows written but not yet checkpointed must be discarded on resume, not
// appended to. This simulates the crash window between the catalog fsync and
// the checkpoint rename.
TEST(Campaign, ResumeTruncatesRowsWrittenPastTheCheckpoint) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "torn.csv");
  o.threads = 2;
  o.stop_after_batches = 3;
  const CampaignStats first = Campaign(o).run();
  ASSERT_FALSE(first.complete);

  // Append junk, as a crashed run that had flushed a batch would have left.
  {
    std::ofstream out(o.out, std::ios::app);
    out << "999,1,2,3,4,5,6,7,0\n";
  }

  CampaignOptions more = o;
  more.resume = true;
  more.stop_after_batches = 0;
  const CampaignStats last = Campaign(more).run();
  EXPECT_TRUE(last.complete);

  for (const std::string &line : ReadLines(o.out))
    EXPECT_EQ(line.rfind("999,", 0), std::string::npos)
        << "the un-checkpointed row survived: " << line;
  EXPECT_EQ(last.solutions, 2u);
}

// Resuming a campaign that is already finished must report it as finished and
// change nothing. The batch loop never runs in that case, so anything the run
// derives only from executing a batch is silently wrong -- which is exactly how
// segments_done was lost the first time this was written.
TEST(Campaign, ResumingACompleteCampaignIsANoOp) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "done.csv");
  o.threads = 4;
  const CampaignStats first = Campaign(o).run();
  ASSERT_TRUE(first.complete);
  const std::vector<std::string> before = ReadLines(o.out);

  CampaignOptions again = o;
  again.resume = true;
  const CampaignStats second = Campaign(again).run();

  EXPECT_TRUE(second.complete) << "an operator would re-run this forever";
  EXPECT_EQ(second.segments_done, first.segments_done);
  EXPECT_EQ(second.norms, first.norms);
  EXPECT_EQ(second.solutions, first.solutions);
  EXPECT_EQ(ReadLines(o.out), before);
}

TEST(Campaign, RefusesToResumeADifferentBound) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "mix.csv");
  o.threads = 2;
  o.stop_after_batches = 1;
  Campaign(o).run();

  CampaignOptions other = o;
  other.bound = 500'000;
  other.resume = true;
  other.stop_after_batches = 0;
  // Merging two bounds into one catalog would produce a file that is neither
  // an exhaustive search to 5e5 nor to 1e6.
  EXPECT_THROW(Campaign(other).run(), std::runtime_error);
}

TEST(Campaign, RejectsANonPositiveCheckpointInterval) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "bad.csv");
  o.checkpoint_seconds = 0;
  EXPECT_THROW(Campaign{o}, std::invalid_argument);
}

// Batch boundaries are a scheduling detail, not part of the answer: rows are
// written in segment order regardless of how segments were chunked. Adaptive
// batching depends on this, since its boundaries vary with machine speed.
TEST(Campaign, CatalogIsIndependentOfBatching) {
  const TempDir dir;
  std::vector<std::string> reference;
  int i = 0;
  for (const std::size_t batch : {std::size_t{1}, std::size_t{7},
                                  std::size_t{4096}, std::size_t{0}}) {
    CampaignOptions o = Base(dir, "b.csv");
    o.out = dir.file(("b" + std::to_string(i++) + ".csv").c_str());
    o.batch = batch;
    o.threads = 4;
    const CampaignStats s = Campaign(o).run();
    EXPECT_EQ(s.norms, 11'083u) << "batch = " << batch;
    const std::vector<std::string> lines = ReadLines(o.out);
    if (reference.empty())
      reference = lines;
    else
      EXPECT_EQ(lines, reference) << "batch = " << batch;
  }
}

// The campaign is the assembled pipeline, so cross-check it against the parts
// driven directly. Different code path, same answer.
TEST(Campaign, AgreesWithTheDirectlyDrivenPipeline) {
  const TempDir dir;
  CampaignOptions o = Base(dir, "x.csv");
  o.threads = 4;
  const CampaignStats s = Campaign(o).run();

  const pte::cpu::Generator g(o.bound);
  pte::cpu::MeetScratch<pte::cpu::trunc_key_t> scratch;
  std::size_t norms = 0, candidates = 0, solutions = 0;
  g.for_each_norm([&](const pte::cpu::NormView &norm) {
    ++norms;
    pte::cpu::meet<pte::cpu::trunc_key_t>(
        norm.D(), norm.nu.data(), scratch, [&](const pte::cpu::E4Sol &hit) {
          ++candidates;
          if (!pte::cpu::confirm_candidate(norm.nu, hit))
            return;
          std::vector<pte::cpu::factor_t> primes = g.primes_of(norm);
          if (!pte::cpu::is_primitive(norm.points, primes, hit))
            return;
          ++solutions;
        });
  });
  EXPECT_EQ(s.norms, norms);
  EXPECT_EQ(s.candidates, candidates);
  EXPECT_EQ(s.solutions, solutions);
}

} // namespace
