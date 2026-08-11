#include "cpu/tower/tower.h"

#include <algorithm>
#include <bit>
#include <cstring>
#include <string>
#include <utility>

namespace pte::cpu {
namespace {

constexpr std::uint32_t kWordBits = 64;

// A sanity bound on the modulus, not a memory bound: the real limit is
// states * m/8 bytes, which is quadratic in m and bites long before this does.
constexpr std::uint32_t kMaxModulus = 1u << 20;

std::size_t words_for(std::uint32_t m) {
  return (m + kWordBits - 1) / kWordBits;
}

bool test_bit(const std::uint64_t *w, std::uint32_t i) {
  return ((w[i >> 6] >> (i & 63)) & 1u) != 0;
}

void set_bit(std::uint64_t *w, std::uint32_t i) {
  w[i >> 6] |= std::uint64_t{1} << (i & 63);
}

bool is_prime(std::uint32_t n) {
  if (n < 2)
    return false;
  for (std::uint64_t d = 2; d * d <= n; ++d)
    if (n % d == 0)
      return false;
  return true;
}

std::uint32_t checked_pow(std::uint32_t p, std::uint32_t k) {
  std::uint64_t m = 1;
  for (std::uint32_t i = 0; i < k; ++i) {
    m *= p;
    if (m > kMaxModulus)
      throw std::invalid_argument("tower: p^k exceeds " +
                                  std::to_string(kMaxModulus));
  }
  return static_cast<std::uint32_t>(m);
}

// The set of distinct states, interned so that the DP's layer is a list of ids
// with multiplicities rather than a map keyed by an m-bit mask. Open addressing
// with a 64-bit fingerprint and a full memcmp on a fingerprint hit -- the same
// shape as gpu/champion's table, and for the same reason: the fingerprint makes
// the common case one compare, and the memcmp makes the answer exact.
class StateTable {
public:
  static constexpr std::uint32_t kMissing = 0xffffffffu;

  explicit StateTable(std::size_t words) : words_(words) { rehash(256); }

  void clear() {
    arena_.clear();
    fingerprints_.clear();
    std::fill(slots_.begin(), slots_.end(), 0u);
    count_ = 0;
  }

  std::uint32_t intern(const std::uint64_t *mask) { return probe(mask, true); }
  std::uint32_t find(const std::uint64_t *mask) const {
    return const_cast<StateTable *>(this)->probe(mask, false);
  }

  std::uint32_t size() const { return count_; }
  const std::uint64_t *mask(std::uint32_t id) const {
    return arena_.data() + static_cast<std::size_t>(id) * words_;
  }

private:
  std::uint64_t hash(const std::uint64_t *mask) const {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (std::size_t w = 0; w < words_; ++w) {
      h ^= mask[w];
      h *= 0x100000001b3ull;
    }
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33;
    return h;
  }

  void rehash(std::size_t slots) {
    slots_.assign(slots, 0u);
    slot_mask_ = slots - 1;
    for (std::uint32_t id = 0; id < count_; ++id) {
      std::uint64_t s = fingerprints_[id] & slot_mask_;
      while (slots_[s] != 0)
        s = (s + 1) & slot_mask_;
      slots_[s] = id + 1;
    }
  }

  std::uint32_t probe(const std::uint64_t *mask, bool insert) {
    const std::uint64_t h = hash(mask);
    std::uint64_t s = h & slot_mask_;
    for (;;) {
      const std::uint32_t entry = slots_[s];
      if (entry == 0) {
        if (!insert)
          return kMissing;
        const std::uint32_t id = count_;
        arena_.insert(arena_.end(), mask, mask + words_);
        fingerprints_.push_back(h);
        slots_[s] = id + 1;
        ++count_;
        if (count_ * 2 >= slots_.size())
          rehash(slots_.size() * 2);
        return id;
      }
      const std::uint32_t id = entry - 1;
      if (fingerprints_[id] == h &&
          std::memcmp(this->mask(id), mask, words_ * sizeof(std::uint64_t)) == 0)
        return id;
      s = (s + 1) & slot_mask_;
    }
  }

  std::size_t words_;
  std::vector<std::uint64_t> arena_;
  std::vector<std::uint64_t> fingerprints_;
  std::vector<std::uint32_t> slots_;
  std::uint32_t count_ = 0;
  std::uint64_t slot_mask_ = 0;
};

// `dbl` holds the mask twice end to end, so every cyclic rotation is a plain
// bit-window extraction rather than a shift-and-or with a wrap-around case.
void build_doubled(const std::uint64_t *mask, std::uint32_t m, std::size_t words,
                   std::uint64_t *dbl, std::size_t dbl_words) {
  std::fill(dbl, dbl + dbl_words, 0ull);
  std::memcpy(dbl, mask, words * sizeof(std::uint64_t));
  const std::size_t q = m >> 6;
  const unsigned r = m & 63;
  for (std::size_t w = 0; w < words; ++w) {
    dbl[q + w] |= mask[w] << r;
    if (r != 0)
      dbl[q + w + 1] |= mask[w] >> (kWordBits - r);
  }
}

// out[i] = mask[(i + shift) mod m], i.e. the set { s - shift : s in mask }.
void rotate_left(const std::uint64_t *dbl, std::uint32_t shift,
                 std::size_t words, std::uint64_t tail, std::uint64_t *out) {
  for (std::size_t w = 0; w < words; ++w) {
    const std::size_t bit = static_cast<std::size_t>(shift) + w * kWordBits;
    const std::size_t q = bit >> 6;
    const unsigned r = bit & 63;
    std::uint64_t v = dbl[q] >> r;
    if (r != 0)
      v |= dbl[q + 1] << (kWordBits - r);
    out[w] = v;
  }
  out[words - 1] &= tail;
}

std::uint64_t popcount(const std::uint64_t *mask, std::size_t words) {
  std::uint64_t n = 0;
  for (std::size_t w = 0; w < words; ++w)
    n += static_cast<std::uint64_t>(std::popcount(mask[w]));
  return n;
}

} // namespace

TowerCounts tower_counts(const TowerParams &params) {
  const std::uint32_t p = params.p;
  const std::uint32_t k = params.k;
  const std::uint32_t depth = params.depth;
  if (!is_prime(p))
    throw std::invalid_argument("tower: p must be prime");
  if (k == 0 || depth == 0)
    throw std::invalid_argument("tower: k and depth must be positive");

  const std::uint32_t m = checked_pow(p, k);
  const std::size_t words = words_for(m);
  const std::uint64_t tail =
      (m % kWordBits) != 0
          ? (std::uint64_t{1} << (m % kWordBits)) - 1
          : ~std::uint64_t{0};
  const std::uint32_t cap = std::min(params.max_valuation, k);

  // The exact tuple count is bounded by m^depth, so refuse the run rather than
  // wrap. Widening the accumulator is a one-line change if the density branch
  // ever needs a modulus this large *and* exact counts; the obstruction branch
  // does not count at all.
  bool counting = params.count_towers;
  if (counting) {
    std::uint64_t bound = 1;
    for (std::uint32_t i = 0; i < depth; ++i) {
      if (bound > (~std::uint64_t{0}) / m)
        throw TowerTooLarge("tower: exact counting would overflow 64 bits at "
                            "p^k = " +
                            std::to_string(m) + ", depth " +
                            std::to_string(depth) +
                            " (set count_towers = false)");
      bound *= m;
    }
  }

  // v_p(t) <= cap. With v_p(0) = k this is "p^{cap+1} does not divide t", and
  // cap == k admits everything including 0.
  std::vector<std::uint8_t> allowed(m, 1u);
  if (cap < k) {
    const std::uint32_t step = checked_pow(p, cap + 1);
    for (std::uint32_t t = 0; t < m; t += step)
      allowed[t] = 0u;
  }
  std::vector<std::uint32_t> squares(m);
  for (std::uint32_t t = 0; t < m; ++t)
    squares[t] = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(t) * t) % m);

  StateTable cur(words), nxt(words);
  std::vector<std::uint64_t> cur_mult, nxt_mult;

  std::vector<std::uint64_t> scratch(words), empty(words, 0ull), rot(words);
  const std::size_t dbl_words = 2 * words + 2;
  std::vector<std::uint64_t> dbl(dbl_words);

  std::vector<std::uint64_t> full(words, ~std::uint64_t{0});
  full[words - 1] &= tail;
  cur.intern(full.data());
  if (counting)
    cur_mult.push_back(1);

  TowerCounts out;
  for (std::uint32_t level = 0; level < depth; ++level) {
    out.states.push_back(cur.size());
    if (params.count_orbits) {
      std::vector<std::uint8_t> seen(cur.size(), 0u);
      std::uint64_t orbits = 0, misses = 0;
      for (std::uint32_t id = 0; id < cur.size(); ++id) {
        if (seen[id] != 0u)
          continue;
        ++orbits;
        build_doubled(cur.mask(id), m, words, dbl.data(), dbl_words);
        for (std::uint32_t shift = 0; shift < m; ++shift) {
          rotate_left(dbl.data(), shift, words, tail, rot.data());
          const std::uint32_t other = cur.find(rot.data());
          if (other == StateTable::kMissing)
            ++misses;
          else
            seen[other] = 1u;
        }
      }
      out.orbits.push_back(orbits);
      out.orbit_misses.push_back(misses);
    }

    nxt.clear();
    nxt_mult.clear();
    std::uint64_t towers = 0, splitting = 0;

    for (std::uint32_t id = 0; id < cur.size(); ++id) {
      const std::uint64_t *tset = cur.mask(id);
      const std::uint64_t mult = counting ? cur_mult[id] : 0;

      // S = { t^2 : t in T, -t in T, v_p(t) <= cap }. A = T n -T n allowed is
      // symmetric and t, -t have the same square, so half the residues suffice.
      std::fill(scratch.begin(), scratch.end(), 0ull);
      for (std::uint32_t t = 0; 2 * t <= m; ++t) {
        if (allowed[t] == 0u || !test_bit(tset, t))
          continue;
        const std::uint32_t neg = (m - t) % m;
        if (!test_bit(tset, neg))
          continue;
        set_bit(scratch.data(), squares[t]);
      }

      // 0 in T_{j+1} iff L_{j+1} in S, so |S| is exactly the number of L that
      // complete this prefix into a splitting tower.
      const std::uint64_t reachable = popcount(scratch.data(), words);
      if (reachable == 0) {
        // Empty S sends all m values of L to the empty set, which is absorbing
        // but still a state -- it is the "+1" in every measured state count.
        const std::uint32_t nid = nxt.intern(empty.data());
        if (counting) {
          if (nid == nxt_mult.size())
            nxt_mult.push_back(0);
          nxt_mult[nid] += mult * m;
        }
        continue;
      }
      ++splitting;
      if (counting)
        towers += mult * reachable;

      build_doubled(scratch.data(), m, words, dbl.data(), dbl_words);
      for (std::uint32_t ell = 0; ell < m; ++ell) {
        rotate_left(dbl.data(), ell, words, tail, rot.data());
        const std::uint32_t nid = nxt.intern(rot.data());
        if (counting) {
          if (nid == nxt_mult.size())
            nxt_mult.push_back(0);
          nxt_mult[nid] += mult;
        }
      }
      if (nxt.size() > params.max_states)
        throw TowerTooLarge("tower: state count exceeded max_states = " +
                            std::to_string(params.max_states) + " at p^k = " +
                            std::to_string(m) + ", level " +
                            std::to_string(level + 1));
    }

    if (counting)
      out.towers.push_back(towers);
    out.splitting_states.push_back(splitting);
    std::swap(cur, nxt);
    std::swap(cur_mult, nxt_mult);
  }
  return out;
}

ValuationFloor valuation_floor(std::uint32_t p, std::uint32_t k,
                               std::uint32_t depth, const TowerParams &base) {
  ValuationFloor out;
  out.p = p;
  out.k = k;
  out.depth = depth;
  TowerParams params = base;
  params.p = p;
  params.k = k;
  params.depth = depth;

  // One run at cap = k-1 already answers "is anything nondegenerate visible at
  // this precision?", and at large k that answer is usually no -- so probing it
  // first replaces k ascending runs with one. Only taken when the caller is not
  // collecting towers_by_cap, which needs every cap anyway.
  if (!params.count_towers && k >= 1) {
    params.max_valuation = k - 1;
    const TowerCounts probe = tower_counts(params);
    if (!probe.splits()) {
      out.floor = k;
      out.states = probe.states;
      return out;
    }
  }

  for (std::uint32_t cap = 0;; ++cap) {
    params.max_valuation = cap;
    const TowerCounts counts = tower_counts(params);
    if (params.count_towers)
      out.towers_by_cap.push_back(counts.towers.back());
    if (counts.splits() || cap == k) {
      out.floor = cap;
      out.states = counts.states;
      return out;
    }
  }
}

} // namespace pte::cpu
