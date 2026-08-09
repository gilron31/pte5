#include "cpu/factorbase/factorbase.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace pte::cpu {
namespace {

// On-disk record: see data/README.md. The prime is not stored -- it is
// a^2 + b^2 -- which halves the artifact.
struct Record {
  std::uint32_t a;
  std::uint32_t b;
};
static_assert(sizeof(Record) == 8, "artifact record must be exactly 8 bytes");

// How many split primes to expect below x, for reserve(). pi(x)/2, padded,
// because the x/ln(x) approximation runs low in the range that matters.
std::size_t EstimateSplitPrimes(factor_t p_max) {
  if (p_max < 5)
    return 0;
  const double x = static_cast<double>(p_max);
  return static_cast<std::size_t>(1.15 * x / std::log(x) / 2.0);
}

// RAII wrapper so a throw mid-write does not leak the handle.
class File {
public:
  File(const std::string &path, const char *mode) : f_(std::fopen(path.c_str(), mode)) {
    if (f_ == nullptr)
      throw std::runtime_error("factorbase: cannot open " + path);
  }
  ~File() {
    if (f_ != nullptr)
      std::fclose(f_);
  }
  File(const File &) = delete;
  File &operator=(const File &) = delete;
  std::FILE *get() const { return f_; }

private:
  std::FILE *f_ = nullptr;
};

} // namespace

SegmentedSieve::SegmentedSieve(factor_t p_max, std::size_t segment)
    : p_max_(p_max), segment_(segment) {
  if (segment_ == 0)
    throw std::invalid_argument("SegmentedSieve: segment width must be positive");
  if (p_max_ < 5)
    return;

  // Base primes to sqrt(p_max). This is the only allocation that grows with the
  // bound, and it grows as sqrt: 11,656 primes at P_max = 1.54e10.
  const auto root = static_cast<factor_t>(std::sqrt(static_cast<double>(p_max_))) + 1;
  std::vector<bool> small(root + 1, true);
  small[0] = false;
  if (root >= 1)
    small[1] = false;
  for (factor_t i = 2; i * i <= root; ++i) {
    if (small[i]) {
      for (factor_t j = i * i; j <= root; j += i)
        small[j] = false;
    }
  }
  for (factor_t i = 2; i <= root; ++i) {
    if (small[i])
      base_.push_back(static_cast<std::uint32_t>(i));
  }
}

FactorBase FactorBase::Build(factor_t p_max) {
  FactorBase fb;
  fb.p_max_ = p_max;
  const std::size_t expect = EstimateSplitPrimes(p_max);
  fb.primes_.reserve(expect);
  fb.points_.reserve(expect);

  const SegmentedSieve sieve(p_max);
  sieve.for_each_split_prime([&](factor_t p) {
    fb.primes_.push_back(p);
    fb.points_.push_back(decompose_prime(p)); // throws if p is not split
  });
  return fb;
}

std::size_t FactorBase::WriteArtifact(const std::string &path, factor_t p_max) {
  const File out(path, "wb");

  // Stream in blocks so the base never has to be resident: this is what lets a
  // 2.6 GB artifact be produced in 0.2 MB of working memory.
  std::vector<Record> buffer;
  buffer.reserve(1u << 16);
  std::size_t written = 0;

  const auto flush = [&]() {
    if (buffer.empty())
      return;
    if (std::fwrite(buffer.data(), sizeof(Record), buffer.size(), out.get()) !=
        buffer.size())
      throw std::runtime_error("factorbase: short write to " + path);
    written += buffer.size();
    buffer.clear();
  };

  const SegmentedSieve sieve(p_max);
  sieve.for_each_split_prime([&](factor_t p) {
    const Gauss z = decompose_prime(p);
    buffer.push_back({static_cast<std::uint32_t>(z.re),
                      static_cast<std::uint32_t>(z.im)});
    if (buffer.size() == buffer.capacity())
      flush();
  });
  flush();
  return written;
}

FactorBase FactorBase::Load(const std::string &path, factor_t p_max) {
  const File in(path, "rb");
  if (std::fseek(in.get(), 0, SEEK_END) != 0)
    throw std::runtime_error("factorbase: cannot seek " + path);
  const long bytes = std::ftell(in.get());
  if (bytes < 0)
    throw std::runtime_error("factorbase: cannot size " + path);
  if (static_cast<std::size_t>(bytes) % sizeof(Record) != 0) {
    throw std::runtime_error("factorbase: " + path + " is " +
                             std::to_string(bytes) +
                             " bytes, not a multiple of the 8-byte record");
  }
  std::rewind(in.get());

  const auto count = static_cast<std::size_t>(bytes) / sizeof(Record);
  FactorBase fb;
  fb.primes_.reserve(count);
  fb.points_.reserve(count);

  std::vector<Record> buffer(1u << 16);
  std::size_t remaining = count;
  factor_t previous = 0;
  bool truncated = false;

  while (remaining > 0 && !truncated) {
    const std::size_t want = std::min(remaining, buffer.size());
    if (std::fread(buffer.data(), sizeof(Record), want, in.get()) != want)
      throw std::runtime_error("factorbase: short read from " + path);
    remaining -= want;

    for (std::size_t i = 0; i < want; ++i) {
      const auto a = static_cast<coord_t>(buffer[i].a);
      const auto b = static_cast<coord_t>(buffer[i].b);
      const auto p = static_cast<factor_t>(a) * a + static_cast<factor_t>(b) * b;

      // Structural validation. A corrupt or misaligned artifact would otherwise
      // feed silently wrong primes into the enumeration, which is exactly the
      // failure an exclusion claim cannot survive.
      if (a < b || b <= 0) {
        throw std::runtime_error("factorbase: " + path + " record " +
                                 std::to_string(count - remaining - want + i) +
                                 " violates a >= b > 0");
      }
      if (p <= previous) {
        throw std::runtime_error("factorbase: " + path +
                                 " is not ascending at p = " +
                                 std::to_string(p));
      }
      if (p_max != 0 && p > p_max) {
        truncated = true;
        break;
      }
      previous = p;
      fb.primes_.push_back(p);
      fb.points_.push_back(Gauss{a, b});
    }
  }

  if (p_max != 0 && !truncated && !fb.primes_.empty() &&
      fb.primes_.back() < p_max) {
    throw std::runtime_error(
        "factorbase: " + path + " reaches only " +
        std::to_string(fb.primes_.back()) + ", short of the requested bound " +
        std::to_string(p_max));
  }
  fb.p_max_ = p_max != 0 ? p_max : (fb.primes_.empty() ? 0 : fb.primes_.back());
  return fb;
}

} // namespace pte::cpu
