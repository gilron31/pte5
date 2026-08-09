// The candidate quadruple, and the canonical form that makes two of them
// comparable.
//
// This lives in core rather than next to any one producer because four
// independent components speak it: the oracle (cpu/enumerate), the production
// meet (cpu/meet), the confirmation step, and the primitivity rule. Putting it
// beside the oracle would force the production meet to include the oracle's
// header -- which drags in the whole Enumerator, and eventually the device
// build, for a 4-integer struct.
#pragma once

#include <tuple>
#include <utility>

namespace pte::cpu {

// Four indices into a point set, read as {a,b} | {c,d}.
struct E4Sol {
  int a = 0, b = 0, c = 0, d = 0;

  // Canonical form: each pair sorted, then the two pairs ordered. A meet may
  // report the same solution under either pair ordering -- the oracle walks
  // 4-subsets, the hashed path walks a sorted key list -- so canonizing is what
  // lets two implementations be compared, sorted and deduplicated at all.
  //
  // It deliberately does *not* merge different pairings of the same 4-subset:
  // {a,b}|{c,d} and {a,c}|{b,d} are genuinely different solutions, and a
  // 4-subset can satisfy both.
  E4Sol canonical() const {
    int p = a, q = b, r = c, s = d;
    if (p > q)
      std::swap(p, q);
    if (r > s)
      std::swap(r, s);
    if (p > r || (p == r && q > s)) {
      std::swap(p, r);
      std::swap(q, s);
    }
    return {p, q, r, s};
  }

  friend bool operator==(const E4Sol &x, const E4Sol &y) {
    return x.a == y.a && x.b == y.b && x.c == y.c && x.d == y.d;
  }
  friend bool operator<(const E4Sol &x, const E4Sol &y) {
    return std::tie(x.a, x.b, x.c, x.d) < std::tie(y.a, y.b, y.c, y.d);
  }
};

} // namespace pte::cpu
