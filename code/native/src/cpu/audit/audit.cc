#include "cpu/audit/audit.h"

#include "cpu/confirm/confirm.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/meet/meet.h"

#include <algorithm>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace pte::cpu {
namespace {

// Trial division. Only ever applied to a handful of known solutions, so the
// naive form is right: anything cleverer would be code that has to be trusted.
std::vector<factor_t> FactorIntoSplitPrimes(norm_t n) {
  std::vector<factor_t> out;
  auto v = static_cast<factor_t>(n);
  for (factor_t d = 2; d * d <= v; ++d) {
    while (v % d == 0) {
      out.push_back(d);
      v /= d;
    }
  }
  if (v > 1)
    out.push_back(v);
  return out;
}

std::vector<E4Sol> Canonical(std::vector<E4Sol> v) {
  for (E4Sol &s : v)
    s = s.canonical();
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

// splitmix64: a deterministic, well-distributed sampler with no dependency on
// the standard library's implementation-defined engines, so a sampled audit is
// reproducible across compilers as well as across runs.
std::uint64_t Mix(std::uint64_t &state) {
  state += 0x9e3779b97f4a7c15ull;
  std::uint64_t z = state;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
}

} // namespace

std::vector<KnownSolution> known_solutions() {
  return {
      {67'405, "bremner F2 k=5 / gilro sporadic"},
      {600'445, "bremner F1 k=1 / gilro sporadic"},
      {5'915'065, "bremner F2 k=6 / gilro sporadic"},
      {55'445'869, "bremner sporadic"},
      {81'524'605, "bremner sporadic / gilro sporadic"},
      {166'045'477, "bremner sporadic / gilro sporadic"},
      {242'048'509, "bremner sporadic"},
      {532'946'005, "bremner sporadic / gilro sporadic"},
      {33'065'996'717LL, "gilro sporadic (paper.md 7 exemplar)"},
      {43'937'864'165LL, "bremner F2 k=7"},
      // F2 k=8 = 2,703,273,548,105 sits above the design bound and is left out
      // deliberately: it is the negative control, and an exhaustive run to
      // 10^12 that reported it would be reporting something out of range.
  };
}

KnownCheck verify_known_solutions(norm_t bound) {
  KnownCheck out;
  MeetScratch<trunc_key_t> scratch;

  for (const KnownSolution &known : known_solutions()) {
    if (known.l1 > bound)
      continue;
    ++out.in_range;

    // Rebuild the norm from scratch -- factor it, decompose each prime, and
    // construct the point set -- so this path shares nothing with the campaign's
    // factor base or its spine.
    const std::vector<factor_t> primes = FactorIntoSplitPrimes(known.l1);
    std::vector<Gauss> points;
    bool split = true;
    for (const factor_t p : primes) {
      if (p % 4 != 1) {
        split = false;
        break;
      }
    }
    if (!split) {
      out.missing.push_back(known.l1);
      continue;
    }
    points = point_set(primes);
    const std::vector<norm_t> nu = nu_values(points);

    bool certified = false;
    meet<trunc_key_t>(static_cast<int>(nu.size()), nu.data(), scratch,
                      [&](const E4Sol &hit) {
                        if (certified)
                          return;
                        if (!confirm_candidate(nu, hit))
                          return;
                        if (!is_primitive(points, primes, hit))
                          return;
                        certified = certify(quadruple(points, hit)).valid();
                      });
    if (certified)
      ++out.found;
    else
      out.missing.push_back(known.l1);
  }
  return out;
}

SegmentCheck verify_segments(const Generator &g, std::size_t count,
                             std::uint64_t seed) {
  SegmentCheck out;
  const std::size_t segments = g.segment_count();
  if (segments == 0 || count == 0)
    return out;

  std::uint64_t state = seed;
  std::set<std::size_t> sample;
  // Sample without replacement, but bail out rather than spin if `count`
  // approaches the population.
  for (std::size_t tries = 0; sample.size() < std::min(count, segments) &&
                              tries < 20 * count + 64;
       ++tries) {
    sample.insert(static_cast<std::size_t>(Mix(state) % segments));
  }

  Walker walker(g);
  MeetScratch<trunc_key_t> scratch;

  for (const std::size_t s : sample) {
    ++out.segments;
    walker.walk(s, [&](const NormView &norm) {
      ++out.norms;
      const std::vector<factor_t> primes = g.primes_of(norm);

      // The reference point set: decompose_prime run again per factor, then the
      // multiplicative rebuild. Independent of the loaded or mapped artifact,
      // so this audits the factor base too.
      const std::vector<Gauss> reference = point_set(primes);
      const std::vector<Gauss> produced(norm.points.begin(), norm.points.end());
      if (reference != produced) {
        out.disagreements.push_back("L1 = " + std::to_string(norm.l1) +
                                    ": point sets differ (" +
                                    std::to_string(produced.size()) + " vs " +
                                    std::to_string(reference.size()) + ")");
        return;
      }

      const std::vector<norm_t> ref_nu = nu_values(reference);
      std::vector<E4Sol> ref_sols;
      for (const E4Sol &q : meet_oracle(ref_nu)) {
        if (is_primitive(reference, primes, q))
          ref_sols.push_back(q);
      }

      std::vector<E4Sol> got_sols;
      meet<trunc_key_t>(norm.D(), norm.nu.data(), scratch,
                        [&](const E4Sol &hit) {
                          if (!confirm_candidate(norm.nu, hit))
                            return;
                          if (!is_primitive(norm.points, primes, hit))
                            return;
                          got_sols.push_back(hit);
                        });

      const std::vector<E4Sol> a = Canonical(std::move(ref_sols));
      const std::vector<E4Sol> b = Canonical(std::move(got_sols));
      if (a != b) {
        out.disagreements.push_back(
            "L1 = " + std::to_string(norm.l1) + ": reference found " +
            std::to_string(a.size()) + " solutions, production found " +
            std::to_string(b.size()));
        return;
      }
      out.solutions += b.size();
    });
  }
  return out;
}

} // namespace pte::cpu
