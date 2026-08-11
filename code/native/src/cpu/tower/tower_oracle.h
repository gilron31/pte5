// The tree condition itself, evaluated one L-tuple at a time.
//
// tower.h computes the whole of (Z/p^k)^d at once, by a forward recursion on
// *sets* with interning and rotation. This file does the opposite in every
// respect: it takes one tuple and walks paper.md 2's binary tree of nested
// radicals downward, enumerating square roots at each node. No sets, no
// interning, no rotations, no multiplicities -- and no memoization either, so
// what it evaluates is the definition and nothing else.
//
// That is what makes it usable as an oracle. tower_count_brute() runs it over
// every tuple in (Z/p^k)^depth and must agree with tower_counts().towers.back()
// exactly, which is the fixture drafts/gpu_implementation_plan.md 6 asks for.
//
// It also answers a question the set DP structurally cannot: **what does one
// specific tuple do?** That is the T4 regression. The section-7 exemplar is a
// real E_4, so its reduction must appear as a tower mod p^k once k exceeds the
// valuations a real solution needs -- and `tuple_valuation_floor` recovers
// those valuations from the reduced tuple alone, with no big-integer arithmetic
// anywhere. Measured, it returns min(ht_p, k): 6 at p = 3 (for k >= 7), 2 at
// p = 7 (for k >= 3), matching the node valuations computed independently in
// drafts/probe_local_model.py. A depth-4 null from the DP at a precision below
// those numbers is therefore correct rather than a bug -- which is the entire
// reason drafts/local_model_plan.md 7 refuses to report one.
#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace pte::cpu {

class TowerOracle {
public:
  // Builds the square-root table for Z/p^k, so the cost is amortised over
  // however many tuples are then asked about.
  TowerOracle(std::uint32_t p, std::uint32_t k);

  std::uint32_t modulus() const { return m_; }

  // Does (L_1..L_d) admit a depth-d tower over Z/p^k with every node value of
  // valuation <= max_valuation? Entries of `ells` are reduced mod p^k first.
  bool splits(std::span<const std::uint32_t> ells,
              std::uint32_t max_valuation) const;

  // The least cap for which it does; k when only degenerate towers exist.
  std::uint32_t tuple_valuation_floor(std::span<const std::uint32_t> ells) const;

private:
  std::uint32_t p_, k_, m_;
  std::vector<std::uint32_t> root_offset_; // CSR over targets in Z/m
  std::vector<std::uint32_t> roots_;
  std::vector<std::uint8_t> valuation_; // v_p, capped at k
};

// Exhaustive: every tuple in (Z/p^k)^depth, evaluated by TowerOracle. m^depth
// tuples, so this is a test fixture and nothing else.
std::uint64_t tower_count_brute(std::uint32_t p, std::uint32_t k,
                                std::uint32_t depth,
                                std::uint32_t max_valuation);

// The section-7 exemplar E_4 = (L_1..L_4). Decimal strings because L_4 has 82
// digits: reducing a decimal string mod m is four lines of Horner, and pulling
// GMP into this directory to hold a constant would not be.
extern const char *const kExemplarE4[4];

// Horner over the digits. Throws std::invalid_argument on a non-digit.
std::uint32_t reduce_decimal(std::string_view digits, std::uint32_t m);

// kExemplarE4 reduced mod m.
std::vector<std::uint32_t> exemplar_e4_mod(std::uint32_t m);

} // namespace pte::cpu
