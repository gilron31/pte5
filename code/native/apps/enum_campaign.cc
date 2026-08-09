// The exhaustive E_4 campaign: enumerate every admissible norm below a bound,
// meet, confirm, certify, and write a catalog.
//
// This is the exhaustive-to-B regime of paper.md 6.2 -- the one that licenses
// an exclusion claim -- so it is checkpointed and resumable by construction.
// Kill it at any point and re-run with --resume; the catalog is identical to
// what an uninterrupted run would have produced.
//
//   ./enum_campaign --bound=10000000000 --factor-base=../../data/factorbase_1e12.bin
//   ./enum_campaign --bound=10000000000 --resume
//
// Pure CPU: no CUDA, no GPU. The device's value in this project is the S-D and
// S-E searches, not this one -- see drafts/implementation_design.md 3.1.
#include "cpu/campaign/campaign.h"

#include "cpu/audit/audit.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <string_view>

namespace {

void Usage() {
  std::printf(
      "usage: enum_campaign [options]\n"
      "\n"
      "  --bound=N            enumerate admissible L1 <= N (required)\n"
      "  --out=PATH           catalog CSV (default catalog.csv)\n"
      "  --factor-base=PATH   load a factor-base artifact instead of sieving\n"
      "  --checkpoint=PATH    checkpoint file (default <out>.ckpt)\n"
      "  --threads=N          worker threads (default: one per core)\n"
      "  --batch=N            fixed segments per checkpoint (default: adaptive)\n"
      "  --checkpoint-secs=N  target seconds between checkpoints (default 5)\n"
      "  --stop-after=N       stop cleanly after N batches, leaving a resumable\n"
      "                       checkpoint (default 0 = run to completion)\n"
      "  --load               read the factor base into the heap instead of\n"
      "                       mmapping it (mmap is the default)\n"
      "  --audit=N            after the run, re-verify N sampled segments\n"
      "                       against the reference implementation\n"
      "  --audit-seed=N       sampling seed for --audit (default 1)\n"
      "  --resume             continue from the checkpoint\n"
      "  --quiet              no progress line\n");
}

// Returns true and fills `value` if `arg` is exactly "--<name>=<value>".
bool Flag(std::string_view arg, std::string_view name, std::string &value) {
  const std::string prefix = "--" + std::string(name) + "=";
  if (arg.substr(0, prefix.size()) != prefix)
    return false;
  value = std::string(arg.substr(prefix.size()));
  return true;
}

double Rate(double count, double seconds) {
  return seconds > 0 ? count / seconds : 0.0;
}

} // namespace

int main(int argc, char **argv) {
  pte::cpu::CampaignOptions opts;
  std::string value;
  std::uint64_t audit_seed = 1;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      Usage();
      return 0;
    }
    if (arg == "--resume") {
      opts.resume = true;
    } else if (arg == "--load") {
      opts.map_factor_base = false;
    } else if (arg == "--quiet") {
      opts.quiet = true;
    } else if (Flag(arg, "bound", value)) {
      opts.bound = std::strtoll(value.c_str(), nullptr, 10);
    } else if (Flag(arg, "out", value)) {
      opts.out = value;
    } else if (Flag(arg, "factor-base", value)) {
      opts.factor_base = value;
    } else if (Flag(arg, "checkpoint", value)) {
      opts.checkpoint = value;
    } else if (Flag(arg, "threads", value)) {
      opts.threads = std::atoi(value.c_str());
    } else if (Flag(arg, "batch", value)) {
      opts.batch = std::strtoull(value.c_str(), nullptr, 10);
    } else if (Flag(arg, "checkpoint-secs", value)) {
      opts.checkpoint_seconds = std::strtod(value.c_str(), nullptr);
    } else if (Flag(arg, "audit", value)) {
      opts.audit_segments = std::strtoull(value.c_str(), nullptr, 10);
    } else if (Flag(arg, "audit-seed", value)) {
      audit_seed = std::strtoull(value.c_str(), nullptr, 10);
    } else if (Flag(arg, "stop-after", value)) {
      opts.stop_after_batches = std::strtoull(value.c_str(), nullptr, 10);
    } else {
      std::fprintf(stderr, "unknown argument: %s\n\n", argv[i]);
      Usage();
      return 2;
    }
  }

  if (opts.bound <= 0) {
    std::fprintf(stderr, "--bound is required and must be positive\n\n");
    Usage();
    return 2;
  }

  try {
    std::printf("E_4 exhaustive campaign\n");
    std::printf("  bound        %lld\n", static_cast<long long>(opts.bound));
    std::printf("  catalog      %s\n", opts.out.c_str());
    std::printf("  factor base  %s%s\n",
                opts.factor_base.empty() ? "(sieved)"
                                         : opts.factor_base.c_str(),
                opts.factor_base.empty()      ? ""
                : opts.map_factor_base        ? " (mmap)"
                                              : " (heap)");

    pte::cpu::Campaign campaign(std::move(opts));
    const pte::cpu::CampaignOptions &o = campaign.options();
    std::printf("  threads      %d\n", o.threads);
    std::printf("  checkpoint   %s%s\n", o.checkpoint.c_str(),
                o.resume ? " (resuming)" : "");
    std::printf("  segments     %zu\n", campaign.generator().segment_count());

    const pte::cpu::CampaignStats s = campaign.run();

    std::printf("\n");
    std::printf("  factor base       %10.2f s\n", s.base_seconds);
    std::printf("  enumeration       %10.2f s\n", s.walk_seconds);
    std::printf("\n");
    std::printf("  segments done     %10zu\n", s.segments_done);
    std::printf("  admissible norms  %10zu   (%.2f M/s)\n", s.norms,
                Rate(static_cast<double>(s.norms), s.walk_seconds) / 1e6);
    std::printf("  meet candidates   %10zu\n", s.candidates);
    std::printf("    rejected, key collision mod 2^64  %10zu\n", s.filtered);
    std::printf("    rejected, imprimitive             %10zu\n",
                s.imprimitive);
    std::printf("  catalogued        %10zu\n", s.solutions);
    std::printf("  UNCERTIFIED       %10zu\n", s.uncertified);
    std::printf("\n");

    // Known solutions are checked on every run: it costs one point-set
    // construction each, and it is the only evidence a run at a bound nobody can
    // exhaust has about its own large-norm behaviour.
    const pte::cpu::KnownCheck known =
        pte::cpu::verify_known_solutions(o.bound);
    std::printf("  known solutions in range %zu, recovered %zu\n",
                known.in_range, known.found);
    for (const auto l1 : known.missing)
      std::printf("    MISSING: L1 = %lld\n", static_cast<long long>(l1));

    bool audit_ok = true;
    if (o.audit_segments > 0) {
      const pte::cpu::SegmentCheck seg = pte::cpu::verify_segments(
          campaign.generator(), o.audit_segments, audit_seed);
      std::printf("  audit: %zu segments, %zu norms, %zu solutions, seed %llu\n",
                  seg.segments, seg.norms, seg.solutions,
                  static_cast<unsigned long long>(audit_seed));
      for (const std::string &d : seg.disagreements)
        std::printf("    DISAGREEMENT: %s\n", d.c_str());
      audit_ok = seg.ok();
    }
    std::printf("\n");

    if (!s.complete) {
      std::printf("  INCOMPLETE -- re-run with --resume to continue.\n");
      return 0;
    }
    if (!known.ok() || !audit_ok) {
      std::printf("  FAILED: the audit did not pass.\n");
      return 1;
    }
    if (s.uncertified > 0) {
      std::printf("  FAILED: %zu catalogued solutions carry no proof.\n",
                  s.uncertified);
      return 1;
    }
    std::printf("  Complete. Every catalogued solution is certified.\n");
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "\ncampaign failed: %s\n", e.what());
    return 1;
  }
}
