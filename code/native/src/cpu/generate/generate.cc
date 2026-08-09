#include "cpu/generate/generate.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace pte::cpu {
namespace {

// Worst-case packing: every level at the point budget, plus the child's
// pre-deduplication scratch, which is twice its parent. Two spare levels cover
// that scratch, so the bound holds without reasoning about how the level sizes
// actually grow.
constexpr std::size_t kBufferPoints =
    static_cast<std::size_t>(kMaxSpineDepth + 2) * kMaxPointCount;

} // namespace

Spine::Spine(const FactorBase &base) : base_(&base) {
  buffer_.resize(kBufferPoints);
  // Level 0 is the empty product: U(1) = {1}, the multiplicative identity the
  // whole construction starts from.
  buffer_[0] = Gauss{1, 0};
  size_[0] = 1;
  offset_[0] = 0;
}

void Spine::push(factor_index_t i) {
  if (depth_ >= kMaxSpineDepth) {
    throw std::length_error(
        "Spine: depth budget " + std::to_string(kMaxSpineDepth) +
        " exceeded; the smallest split prime is 5, so a bound of B needs "
        "log_5(B) levels");
  }

  const auto level = static_cast<std::size_t>(depth_);
  const int parent_size = size_[level];
  const std::size_t start =
      static_cast<std::size_t>(offset_[level]) +
      static_cast<std::size_t>(parent_size);
  if (start + static_cast<std::size_t>(2 * parent_size) > buffer_.size()) {
    throw std::length_error("Spine: point buffer exhausted at depth " +
                            std::to_string(depth_));
  }

  // U(cq) = { canonize(z*q), canonize(z*conj(q)) : z in U(c) }, both products
  // from four multiplies rather than eight.
  const Gauss q = base_->point(i);
  Gauss *const dst = buffer_.data() + start;
  const Gauss *const src = buffer_.data() + offset_[level];
  int m = 0;
  for (int k = 0; k < parent_size; ++k) {
    Gauss a;
    Gauss b;
    src[k].mul_and_conj_mul(q, a, b);
    dst[m++] = a.canonize();
    dst[m++] = b.canonize();
  }

  // The collapse is genuinely many-to-one -- see the header for the three-way
  // case at 125 -- so this is a full sort and unique, not a special case.
  std::sort(dst, dst + m);
  m = static_cast<int>(std::unique(dst, dst + m) - dst);

  if (m > kMaxPointCount) {
    throw std::length_error("Spine: point budget " +
                            std::to_string(kMaxPointCount) + " exceeded (" +
                            std::to_string(m) + " points)");
  }

  path_[level] = i;
  offset_[level + 1] = static_cast<int>(start);
  size_[level + 1] = m;
  value_ *= static_cast<norm_t>(base_->prime(i));
  ++depth_;
}

Generator::Generator(norm_t bound) : bound_(bound) {
  if (bound < 0)
    throw std::invalid_argument("Generator: negative bound " +
                                std::to_string(bound));
  // A bound below the smallest admissible norm is a legitimate query with an
  // empty answer: 5 * 13 * 17 = 1105 is the smallest product of split primes
  // with four points, so "nothing below 1000" is simply true.
  if (bound < 5 * 13 * 17)
    return;
  base_ = FactorBase::Build(static_cast<factor_t>(bound / kSmallestCofactor));
  usable_ = base_.size();
}

Generator::Generator(norm_t bound, FactorBase base)
    : bound_(bound), base_(std::move(base)) {
  if (bound < 0)
    throw std::invalid_argument("Generator: negative bound " +
                                std::to_string(bound));
  if (bound < 5 * 13 * 17) {
    base_ = FactorBase{};
    return;
  }
  const auto needed = static_cast<factor_t>(bound / kSmallestCofactor);
  if (base_.p_max() < needed) {
    throw std::invalid_argument(
        "Generator: factor base reaches only " +
        std::to_string(base_.p_max()) + ", short of the " +
        std::to_string(needed) + " a bound of " + std::to_string(bound) +
        " requires; a short base silently deletes every norm whose largest "
        "prime factor lies above it");
  }
  // The base may reach far past what this bound can use -- the shipped 10^12
  // artifact against a 10^10 campaign, for instance. Find the cut once.
  const auto &primes = base_.primes();
  usable_ = static_cast<std::size_t>(
      std::upper_bound(primes.begin(), primes.end(), needed) - primes.begin());
}

std::vector<factor_t> Generator::primes_of(const NormView &norm) const {
  std::vector<factor_t> out;
  out.reserve(norm.factors.size());
  for (const factor_index_t i : norm.factors)
    out.push_back(base_.prime(i));
  return out;
}

} // namespace pte::cpu
