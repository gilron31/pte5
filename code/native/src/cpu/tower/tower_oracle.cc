#include "cpu/tower/tower_oracle.h"

#include <stdexcept>
#include <string>

namespace pte::cpu {
namespace {

bool is_prime(std::uint32_t n) {
  if (n < 2)
    return false;
  for (std::uint64_t d = 2; d * d <= n; ++d)
    if (n % d == 0)
      return false;
  return true;
}

std::uint32_t pow_checked(std::uint32_t p, std::uint32_t k) {
  std::uint64_t m = 1;
  for (std::uint32_t i = 0; i < k; ++i) {
    m *= p;
    if (m > (1u << 20))
      throw std::invalid_argument("tower_oracle: p^k too large");
  }
  return static_cast<std::uint32_t>(m);
}

} // namespace

const char *const kExemplarE4[4] = {
    "33065996717",
    "489628056848329146064",
    "175480010455650701584492675662518592000000",
    // 40042900368028062136207226327668992000000^2, i.e. L_4 is a perfect square
    // -- which is the top node value sqrt(L_4) being an integer.
    "1603433869883821914072519205964440449600800710383853000439352318296064000"
    "000000000",
};

std::uint32_t reduce_decimal(std::string_view digits, std::uint32_t m) {
  std::uint64_t r = 0;
  for (const char c : digits) {
    if (c < '0' || c > '9')
      throw std::invalid_argument("reduce_decimal: not a decimal digit");
    r = (r * 10 + static_cast<std::uint64_t>(c - '0')) % m;
  }
  return static_cast<std::uint32_t>(r);
}

std::vector<std::uint32_t> exemplar_e4_mod(std::uint32_t m) {
  std::vector<std::uint32_t> ells;
  ells.reserve(4);
  for (const char *const digits : kExemplarE4)
    ells.push_back(reduce_decimal(digits, m));
  return ells;
}

TowerOracle::TowerOracle(std::uint32_t p, std::uint32_t k) : p_(p), k_(k) {
  if (!is_prime(p))
    throw std::invalid_argument("tower_oracle: p must be prime");
  if (k == 0)
    throw std::invalid_argument("tower_oracle: k must be positive");
  m_ = pow_checked(p, k);

  valuation_.resize(m_);
  valuation_[0] = static_cast<std::uint8_t>(k); // v_p(0) = k, capped
  for (std::uint32_t t = 1; t < m_; ++t) {
    std::uint32_t v = 0, x = t;
    while (x % p == 0) {
      x /= p;
      ++v;
    }
    valuation_[t] = static_cast<std::uint8_t>(v);
  }

  // CSR: for each target x, the list of t with t^2 == x. Counting sort, so the
  // roots come out ascending and the table is built in two linear passes.
  std::vector<std::uint32_t> target(m_);
  root_offset_.assign(static_cast<std::size_t>(m_) + 1, 0u);
  for (std::uint32_t t = 0; t < m_; ++t) {
    target[t] =
        static_cast<std::uint32_t>((static_cast<std::uint64_t>(t) * t) % m_);
    ++root_offset_[target[t] + 1];
  }
  for (std::uint32_t x = 0; x < m_; ++x)
    root_offset_[x + 1] += root_offset_[x];
  roots_.resize(m_);
  std::vector<std::uint32_t> fill(root_offset_.begin(), root_offset_.end() - 1);
  for (std::uint32_t t = 0; t < m_; ++t)
    roots_[fill[target[t]]++] = t;
}

bool TowerOracle::splits(std::span<const std::uint32_t> ells,
                         std::uint32_t max_valuation) const {
  const std::uint32_t depth = static_cast<std::uint32_t>(ells.size());
  if (depth == 0)
    return true;
  const std::uint32_t cap = max_valuation < k_ ? max_valuation : k_;

  // node(c, j): does E_j(x) = c admit a full binary tree of node values? The
  // level-j node values are the t with t^2 = L_j + c, and both subtrees -- the
  // one under +t and the one under -t -- must themselves be full.
  const auto node = [&](auto &&self, std::uint32_t c, std::uint32_t j) -> bool {
    if (j == 0)
      return true;
    const std::uint32_t ell = ells[j - 1] % m_;
    const std::uint32_t x = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(ell) + c) % m_);
    for (std::uint32_t i = root_offset_[x]; i < root_offset_[x + 1]; ++i) {
      const std::uint32_t t = roots_[i];
      if (valuation_[t] > cap)
        continue;
      if (self(self, t, j - 1) && self(self, (m_ - t) % m_, j - 1))
        return true;
    }
    return false;
  };
  return node(node, 0u, depth);
}

std::uint32_t
TowerOracle::tuple_valuation_floor(std::span<const std::uint32_t> ells) const {
  for (std::uint32_t cap = 0; cap < k_; ++cap)
    if (splits(ells, cap))
      return cap;
  return k_;
}

std::uint64_t tower_count_brute(std::uint32_t p, std::uint32_t k,
                                std::uint32_t depth,
                                std::uint32_t max_valuation) {
  const TowerOracle oracle(p, k);
  const std::uint32_t m = oracle.modulus();
  std::vector<std::uint32_t> ells(depth, 0u);
  std::uint64_t count = 0;
  for (;;) {
    if (oracle.splits(ells, max_valuation))
      ++count;
    std::uint32_t i = 0;
    for (; i < depth; ++i) {
      if (++ells[i] < m)
        break;
      ells[i] = 0;
    }
    if (i == depth)
      return count;
  }
}

} // namespace pte::cpu
