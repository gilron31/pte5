// The S-E hit: two disjoint pairs of point indices with equal key.
//
// Why this is in core rather than in cpu/champion or gpu/champion, when
// REVIEW.md D5 argued the opposite for MeetHit: the argument there was about a
// *dependency edge*, not about duplication -- sharing pte::cpu::Quadruple would
// have made pte_gpu_meet link pte_cpu, dragging GMP and Boost into the GPU half
// and weakening the oracle that the cross-check rests on. D5 says explicitly
// that "if a third consumer ever wants one type, the right home is a neutral
// src/core/ that both halves depend on -- not cpu/". S-E is that case: the CPU
// reference, the device path and the demo that compares them all speak this
// type, and core carries no GMP and no CUDA, so no edge is created in either
// direction.
//
// Deliberately *not* pte::cpu::E4Sol: indices here are uint32 because a
// champion point set has D ~ 2^19 entries against E4Sol's `int` over D <= 192,
// and the two searches are sized a hundred thousand times apart.
#pragma once

#include <cstdint>

namespace pte::core {

// Indices into a champion norm's point set, read as {a,b} | {c,d}.
struct PairCollision {
  std::uint32_t a = 0, b = 0, c = 0, d = 0;

  // Canonical form: each pair sorted, then the two pairs ordered. Two
  // implementations of the same search legitimately report a solution under
  // either pair ordering -- the CPU reference scans a sorted key list, the
  // device kernel reports (mine, the one already in the table) -- so canonizing
  // is what makes the two comparable at all.
  //
  // Like E4Sol::canonical it does *not* merge different pairings of the same
  // 4-subset: {a,b}|{c,d} and {a,c}|{b,d} are different solutions, and one
  // 4-subset can satisfy both.
  PairCollision canonical() const {
    std::uint32_t p = a, q = b, r = c, s = d;
    if (p > q) {
      const std::uint32_t t = p;
      p = q;
      q = t;
    }
    if (r > s) {
      const std::uint32_t t = r;
      r = s;
      s = t;
    }
    if (p > r || (p == r && q > s))
      return {r, s, p, q};
    return {p, q, r, s};
  }

  friend bool operator==(const PairCollision &x, const PairCollision &y) {
    return x.a == y.a && x.b == y.b && x.c == y.c && x.d == y.d;
  }
  friend bool operator!=(const PairCollision &x, const PairCollision &y) {
    return !(x == y);
  }
  friend bool operator<(const PairCollision &x, const PairCollision &y) {
    if (x.a != y.a)
      return x.a < y.a;
    if (x.b != y.b)
      return x.b < y.b;
    if (x.c != y.c)
      return x.c < y.c;
    return x.d < y.d;
  }
};

} // namespace pte::core
