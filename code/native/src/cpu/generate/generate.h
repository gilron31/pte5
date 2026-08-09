// Norm generation with point sets carried down the DFS spine.
//
// The enumeration walks a tree of cofactors: a node is a product of split
// primes, a child is that product times one more prime copy. The point sets of
// parent and child are related, and the whole reason this file exists is that
// the previous implementation threw that relation away and rebuilt each set from
// {1, 0} on every norm.
//
// ---------------------------------------------------------------------------
// The relation
// ---------------------------------------------------------------------------
//
//   U(c * q) = { canonize(z*q), canonize(z*conj(q)) : z in U(c) },  deduplicated
//
// This is not obvious and is worth the paragraph, because a plausible-looking
// version of it is wrong. U(c) is not the set of Gaussian integers of norm c --
// it is that set modulo units *and conjugation*, which is what canonize()
// projects onto. The identity above holds because the full set W(c) is closed
// under conjugation, so conj(W(c)*q) = W(c)*conj(q): every class of the child
// already has a representative of the form z*q with z in W(c), and the two
// choices z, conj(z) of class representative give the two children z*q and
// z*conj(q). Multiplication by a unit passes straight through, since
// (iz)q = i(zq).
//
// The deduplication is not optional and is not a single collision. At c = 25,
// q = 2+i, the parent U(25) = {(5,0), (4,3)} produces (10,5), (10,5), (10,5)
// and (11,2) -- a three-way collapse -- giving U(125) = {(10,5), (11,2)}. Any
// scheme that assumes "at most one duplicate, from the self-conjugate point"
// silently produces a point set that is too large, and a point set that is too
// large means nu is no longer injective on it, which is the precondition the
// entire meet rests on. So: canonize, sort, unique, at every level.
//
// ---------------------------------------------------------------------------
// The budgets
// ---------------------------------------------------------------------------
//
// Both are exact, from a search over products of split primes <= B, and both
// agree with the exhaustive enumeration where one exists (D_max = 16 at 10^6
// and 36 at 10^8, matching drafts/implementation_design.md 2):
//
//   B        max prime copies (depth)   max divisors   max |U| = D
//   10^6              8                       32            16
//   10^8             11                       72            36
//   10^10            14                      160            80
//   10^12            17                      384           192
//
// Depth is driven by the *smallest* split prime, not by the number of distinct
// primes: 5^17 = 7.6e11 is an admissible norm below 10^12 with seventeen prime
// copies but only nine points. D is driven by the opposite extreme. An earlier
// draft budgeted depth 13 by extrapolating the distinct-prime chain, which
// 5^17 exceeds outright.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/gaussian/gaussian.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace pte::cpu {

// Compile-time capacity, sized for the design bound kMaxNorm = 10^12 with
// headroom over the exact figures above. Exceeding either is a throw, never a
// truncation: a spine that quietly dropped points would delete solutions.
inline constexpr int kMaxSpineDepth = 20;  // 17 needed at 10^12
inline constexpr int kMaxPointCount = 192; // exact at 10^12

// The largest prime that can divide an admissible norm <= B is B / 65: with
// exponent 1 its cofactor needs prod(b_i + 1) >= 4, and the smallest such
// cofactor is 5 * 13 = 65.
//
// Exactly tight, and checked by brute force in enumerate_test.cc rather than
// trusted: at B = 10^6 the largest prime factor that actually occurs is 15,377
// against a cutoff of 15,384, from 5 * 13 * 15377 = 999,505. The reference
// enumerator keeps its own copy of this constant on purpose -- it is the
// independent implementation, and sharing the one number that decides how far
// the factor base reaches would make the cross-check weaker.
inline constexpr norm_t kSmallestCofactor = 65;

// One admissible norm, as seen by a callback. Every span aliases the walk's
// working state and is invalidated as soon as the callback returns.
struct NormView {
  norm_t l1 = 0;
  std::span<const factor_index_t> factors; // ascending, with multiplicity
  std::span<const Gauss> points;           // canonized, sorted, distinct
  std::span<const norm_t> nu;              // nu(z) = Re(z^2), in point order

  // The meet's dispatch bucket. Free here -- the walk already knows it, so
  // nothing has to compute prod(b_i + 1) a second time.
  int D() const { return static_cast<int>(points.size()); }
};

// The DFS stack, carrying one point set per level.
//
// Storage is a single buffer allocated once, with levels packed end to end.
// push() writes the child immediately after the parent, so the two levels the
// inner loop touches are contiguous; pop() is a decrement. Nothing here
// allocates after construction, which is the point: point-set construction was
// the largest single CPU cost in the enumeration, and a per-norm vector is most
// of it.
class Spine {
public:
  explicit Spine(const FactorBase &base);

  // Extend by one copy of base.prime(i). O(D) multiplies plus the sort.
  // Throws std::length_error if the depth or point budget would be exceeded.
  void push(factor_index_t i);

  // Drop the deepest level. Undefined at depth 0, which the DFS cannot reach.
  void pop() { --depth_; }

  int depth() const { return depth_; }

  // The cofactor formed by the current path, and its point set.
  norm_t value() const { return value_; }
  std::span<const Gauss> points() const {
    return {buffer_.data() + offset_[static_cast<std::size_t>(depth_)],
            static_cast<std::size_t>(size_[static_cast<std::size_t>(depth_)])};
  }
  std::span<const factor_index_t> factors() const {
    return {path_.data(), static_cast<std::size_t>(depth_)};
  }

private:
  const FactorBase *base_;
  int depth_ = 0;
  norm_t value_ = 1;
  std::array<factor_index_t, kMaxSpineDepth> path_{};
  std::array<int, kMaxSpineDepth + 1> size_{};
  std::array<int, kMaxSpineDepth + 1> offset_{};
  std::vector<Gauss> buffer_;
};

// Drives the enumeration for one bound over a factor base.
//
// Supersedes Enumerator for production use. Enumerator stays as the reference:
// it rebuilds every point set from scratch, which is slow and obviously
// correct, and generate_test.cc checks this against it norm for norm. Same
// discipline as meet vs meet_oracle.
class Generator {
public:
  // Sieves and decomposes its own factor base.
  explicit Generator(norm_t bound);

  // Takes a factor base the caller already has -- typically loaded from
  // data/factorbase_*.bin, which skips Cornacchia entirely. Throws
  // std::invalid_argument if the base does not reach bound / 65, since a short
  // base silently deletes every norm with a large prime factor.
  Generator(norm_t bound, FactorBase base);

  norm_t bound() const { return bound_; }
  const FactorBase &base() const { return base_; }

  // The enumeration is partitioned into one segment per usable prime: segment s
  // holds exactly the norms whose *largest* prime factor is base.prime(s).
  // Every admissible norm has one, so the segments are disjoint and cover
  // everything -- which is what makes them independently schedulable and
  // independently resumable.
  //
  // Segmenting by the largest prime rather than the smallest is not a
  // preference, it is the difference between a parallel campaign and a serial
  // one. Measured at B = 10^8 (cost = quad-partitions + points):
  //
  //   segmentation      non-empty segments   largest segment   top 16
  //   by smallest prime         42                89.4%         99.9%
  //   by largest prime      58,381                 2.2%         29.8%
  //
  // Smallest-prime segments are what the natural ascending DFS produces, and
  // they cap the speedup at 1.1x on any number of cores.
  std::size_t segment_count() const { return usable_; }

  // Calls cb(const NormView&) once per admissible L1 <= bound(), in no
  // particular order.
  //
  // Admissible means a product of split primes whose point set has at least
  // four elements, i.e. prod(b_i + 1) >= 8. Repeated primes are included:
  // roughly a fifth of admissible norms are non-squarefree and they can carry
  // primitive solutions.
  //
  // Allocates one Walker per call, not per norm, and touches no shared state,
  // so separate threads may call this on the same Generator. A campaign should
  // drive Walker directly instead, one per thread.
  template <typename F> void for_each_norm(F &&cb) const;

  // The prime values behind a NormView's factor indices, with multiplicity.
  // Only needed by is_primitive; the enumeration itself stays on indices.
  std::vector<factor_t> primes_of(const NormView &norm) const;

private:
  norm_t bound_ = 0;
  FactorBase base_;
  // Primes at or below bound / 65, which is where the factor base stops being
  // useful. Stored rather than recomputed because a base loaded from the 10^12
  // artifact is 60x longer than a 10^10 campaign can use, and the DFS's own
  // "value > bound / p" guard does not cut the top level short -- at value = 1
  // it admits every prime up to the bound itself.
  std::size_t usable_ = 0;
};

// One thread's working state for walking segments.
//
// Exists because the Spine's buffer is ~68 KB: allocating one per segment would
// mean 4.3 million allocations over a B = 10^10 campaign. A Walker is built
// once per thread and reused across every segment that thread is handed.
class Walker {
public:
  explicit Walker(const Generator &g)
      : g_(&g), spine_(g.base()), nu_(kMaxPointCount) {}

  // Calls cb(const NormView&) for every admissible norm whose largest prime
  // factor is base.prime(segment). Returns the spine to depth 0, so the same
  // Walker can be handed the next segment immediately.
  template <typename F> void walk(std::size_t segment, F &&cb);

private:
  const Generator *g_;
  Spine spine_;
  // nu is derived at emission rather than inside Spine::push, because push runs
  // on every node of the cofactor tree and only some of those nodes are
  // admissible norms -- computing it there would pay for the ones nobody asks
  // about. Sized once; the meet reads this buffer directly, which is what keeps
  // the campaign's inner loop free of a per-norm vector.
  std::vector<norm_t> nu_;
};

// ---------------------------------------------------------------------------
// Implementation of the templated members
// ---------------------------------------------------------------------------

template <typename F> void Walker::walk(std::size_t segment, F &&cb) {
  if (segment >= g_->segment_count())
    return;
  const norm_t bound = g_->bound();
  const FactorBase &base = g_->base();
  const auto top = static_cast<norm_t>(
      base.prime(static_cast<factor_index_t>(segment)));
  if (top > bound)
    return;

  const auto emit = [&](norm_t v) {
    const std::span<const Gauss> points = spine_.points();
    for (std::size_t k = 0; k < points.size(); ++k)
      nu_[k] = static_cast<norm_t>(points[k].nu());
    cb(NormView{v, spine_.factors(), points,
                std::span<const norm_t>(nu_.data(), points.size())});
  };

  // The cofactor, over primes strictly below the segment's own. Walked in
  // *ascending* index order even though the segment is chosen by the largest
  // prime: ascending is what makes "value > bound / p" a valid break rather
  // than a skip, which is the difference between O(depth) and O(segment) work
  // per node.
  const auto dfs = [&](auto &&self, std::size_t start, norm_t value,
                       std::size_t divisors) -> void {
    for (std::size_t i = start; i < segment; ++i) {
      const auto p =
          static_cast<norm_t>(base.prime(static_cast<factor_index_t>(i)));
      if (value > bound / p)
        break;

      norm_t v = value * p;
      std::size_t exponent = 1;
      while (true) {
        spine_.push(static_cast<factor_index_t>(i));
        const std::size_t d = divisors * (exponent + 1);
        if (d >= 8) // >= 4 points, the minimum for an E_4
          emit(v);
        self(self, i + 1, v, d);

        if (v > bound / p)
          break;
        v *= p;
        ++exponent;
      }
      for (std::size_t e = 0; e < exponent; ++e)
        spine_.pop();
    }
  };

  // The segment's own prime, at every exponent that fits.
  norm_t v = top;
  std::size_t exponent = 1;
  while (true) {
    spine_.push(static_cast<factor_index_t>(segment));
    const std::size_t d = exponent + 1;
    if (d >= 8) // a pure power, p^7 or higher
      emit(v);
    dfs(dfs, 0, v, d);

    if (v > bound / top)
      break;
    v *= top;
    ++exponent;
  }
  for (std::size_t e = 0; e < exponent; ++e)
    spine_.pop();
}

template <typename F> void Generator::for_each_norm(F &&cb) const {
  if (usable_ == 0)
    return;
  Walker walker(*this);
  for (std::size_t s = 0; s < usable_; ++s)
    walker.walk(s, cb);
}

} // namespace pte::cpu
