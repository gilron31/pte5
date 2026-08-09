// The exhaustive campaign: generate, meet, confirm, certify, catalog.
//
// This is the regime paper.md 6.2 calls *exhaustive-to-B*, the one that
// licenses an exclusion claim. Which means the interesting engineering here is
// not throughput -- it is that a run which stops early must be
// indistinguishable, on resumption, from one that never stopped. A campaign
// that silently skipped a segment would produce a catalog that looks exactly
// like a complete one.
//
// ---------------------------------------------------------------------------
// Segments
// ---------------------------------------------------------------------------
//
// The unit of work, of scheduling, and of checkpointing is a *segment*: all
// norms whose largest prime factor is one particular prime (see
// Generator::segment_count). Segments partition the search space, so:
//
//   - any subset of them can run in parallel with no coordination, and
//   - "segments 0..k are done" is a complete description of progress.
//
// Segments are processed in batches. Within a batch, threads pull segments off
// an atomic counter; at the end of a batch the rows are written back **in
// segment order**, flushed, and only then is the checkpoint replaced. So the
// catalog is byte-identical regardless of thread count, and a kill costs at
// most one batch.
//
// ---------------------------------------------------------------------------
// The commit order, which is the whole of the resume story
// ---------------------------------------------------------------------------
//
//   1. write the batch's rows to the catalog
//   2. fflush + fsync the catalog          <- rows are durable
//   3. write checkpoint.tmp, fsync it
//   4. rename checkpoint.tmp -> checkpoint <- atomic on POSIX
//
// A crash between 2 and 4 leaves rows on disk that the checkpoint does not
// know about. Resume therefore *truncates* the catalog back to the checkpointed
// byte offset before appending. Truncating is what makes replay idempotent;
// appending blindly would duplicate a batch.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/generate/generate.h"

#include <cstddef>
#include <string>

namespace pte::cpu {

struct CampaignOptions {
  norm_t bound = 0;

  // 0 means "one per hardware thread".
  int threads = 0;

  // Path to a factor-base artifact (data/factorbase_*.bin). Empty sieves one,
  // which costs 86% more than loading because of Cornacchia.
  std::string factor_base;

  std::string out = "catalog.csv";

  // Defaults to out + ".ckpt".
  std::string checkpoint;

  // Continue from the checkpoint instead of starting over. Refuses to resume a
  // checkpoint written for a different bound.
  bool resume = false;

  // Segments per checkpoint. **0 means adaptive**, which is the default and
  // almost always what you want; a positive value pins the batch size and is
  // mainly useful for tests that need deterministic boundaries.
  //
  // Fixed batches do not work here, and the reason is worth recording. Segment
  // cost falls steeply with index -- a segment's norms are the ones whose
  // largest prime factor is that segment's prime, and small primes admit far
  // more smooth cofactors. Measured over the B = 10^10 campaign with a fixed
  // batch of 4096 segments:
  //
  //   batch 0 (segments      0..4096)   39.4% of all norms
  //   batch 1 (segments   4096..8192)    9.4%
  //   batch 2                            5.2%
  //   ...      1,056 batches in total, the last few thousandths of a percent
  //
  // So the first checkpoint lands 39% of the way in. At B = 10^10 that is 7
  // seconds; at the design bound it would be a quarter of an hour of
  // unprotected work. Geometric batch growth does not fix it either -- the
  // concentration is inside the first block, not spread across the first few.
  std::size_t batch = 0;

  // Target wall-clock seconds between checkpoints, when `batch` is adaptive.
  // The batch size is rescaled after each batch to converge on this, which
  // needs no model of the cost curve and works at any bound.
  double checkpoint_seconds = 5.0;

  // Stop cleanly after this many batches (0 = run to completion), leaving a
  // valid checkpoint. Both a way to bound a session and the hook the
  // resume-equivalence test uses to simulate a kill.
  std::size_t stop_after_batches = 0;

  bool quiet = false;
};

struct CampaignStats {
  std::size_t segments_done = 0;
  std::size_t norms = 0;
  std::size_t candidates = 0;   // meet hits, before exact confirmation
  std::size_t filtered = 0;     // rejected by exact confirmation (mod 2^64 collisions)
  std::size_t imprimitive = 0;  // rejected as scaled duplicates
  std::size_t uncertified = 0;  // reached the catalog but failed the certificate
  std::size_t solutions = 0;    // rows written
  bool complete = false;        // every segment processed

  double base_seconds = 0;
  double walk_seconds = 0;
};

// Runs one campaign. Not copyable; hold it for the duration of a run.
class Campaign {
public:
  explicit Campaign(CampaignOptions opts);
  Campaign(const Campaign &) = delete;
  Campaign &operator=(const Campaign &) = delete;

  CampaignStats run();

  const Generator &generator() const { return generator_; }
  const CampaignOptions &options() const { return opts_; }

  // The catalog's column header, without the trailing newline. Public because
  // both the writer and the resume path need to agree on it.
  static const char *header();

private:
  CampaignOptions opts_;
  // Declared before generator_ so it is alive when the factor base is built:
  // member initialization follows declaration order, and building the base is
  // the one startup cost worth reporting separately.
  double base_seconds_ = 0;
  Generator generator_;
};

} // namespace pte::cpu
