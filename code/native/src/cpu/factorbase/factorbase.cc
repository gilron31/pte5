#include "cpu/factorbase/factorbase.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace pte::cpu {
namespace {

using Record = FactorBase::Record;

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
  File(const std::string &path, const char *mode)
      : f_(std::fopen(path.c_str(), mode)) {
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

factor_t PrimeOf(const Record &r) {
  return static_cast<factor_t>(r.a) * r.a + static_cast<factor_t>(r.b) * r.b;
}

// Structural validation of a record run. A corrupt or misaligned artifact would
// otherwise feed silently wrong primes into the enumeration, which is exactly
// the failure an exclusion claim cannot survive.
//
// Returns the number of records at or below p_max (all of them when p_max is
// 0), throwing on any violation before that point.
std::size_t ValidatePrefix(const Record *records, std::size_t count,
                           factor_t p_max, const std::string &path) {
  factor_t previous = 0;
  for (std::size_t i = 0; i < count; ++i) {
    const Record &r = records[i];
    if (r.a < r.b || r.b == 0) {
      throw std::runtime_error("factorbase: " + path + " record " +
                               std::to_string(i) + " violates a >= b > 0");
    }
    const factor_t p = PrimeOf(r);
    if (p <= previous) {
      throw std::runtime_error("factorbase: " + path +
                               " is not ascending at p = " +
                               std::to_string(p));
    }
    if (p_max != 0 && p > p_max)
      return i;
    previous = p;
  }
  if (p_max != 0 && previous < p_max) {
    throw std::runtime_error("factorbase: " + path + " reaches only " +
                             std::to_string(previous) +
                             ", short of the requested bound " +
                             std::to_string(p_max));
  }
  return count;
}

// Size of a file in whole records, throwing if it is not a whole number.
std::size_t RecordCount(std::size_t bytes, const std::string &path) {
  if (bytes % sizeof(Record) != 0) {
    throw std::runtime_error("factorbase: " + path + " is " +
                             std::to_string(bytes) +
                             " bytes, not a multiple of the 8-byte record");
  }
  return bytes / sizeof(Record);
}

} // namespace

SegmentedSieve::SegmentedSieve(factor_t p_max, std::size_t segment)
    : p_max_(p_max), segment_(segment) {
  if (segment_ == 0)
    throw std::invalid_argument(
        "SegmentedSieve: segment width must be positive");
  if (p_max_ < 5)
    return;

  // Base primes to sqrt(p_max). This is the only allocation that grows with the
  // bound, and it grows as sqrt: 11,656 primes at P_max = 1.54e10.
  const auto root =
      static_cast<factor_t>(std::sqrt(static_cast<double>(p_max_))) + 1;
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

FactorBase::~FactorBase() { Reset(); }

void FactorBase::Reset() {
  if (mapping_ != nullptr) {
    ::munmap(mapping_, mapping_bytes_);
    mapping_ = nullptr;
    mapping_bytes_ = 0;
  }
  records_ = nullptr;
  count_ = 0;
}

FactorBase::FactorBase(FactorBase &&other) noexcept
    : count_(other.count_), p_max_(other.p_max_),
      owned_(std::move(other.owned_)), mapping_(other.mapping_),
      mapping_bytes_(other.mapping_bytes_) {
  // owned_ moved, so its data() is the new home; a mapped base keeps its
  // pointer. Recomputing rather than copying records_ is what makes the move
  // correct for both.
  records_ = mapping_ != nullptr
                 ? static_cast<const Record *>(mapping_)
                 : owned_.data();
  other.mapping_ = nullptr;
  other.mapping_bytes_ = 0;
  other.records_ = nullptr;
  other.count_ = 0;
}

FactorBase &FactorBase::operator=(FactorBase &&other) noexcept {
  if (this != &other) {
    Reset();
    count_ = other.count_;
    p_max_ = other.p_max_;
    owned_ = std::move(other.owned_);
    mapping_ = other.mapping_;
    mapping_bytes_ = other.mapping_bytes_;
    records_ = mapping_ != nullptr ? static_cast<const Record *>(mapping_)
                                   : owned_.data();
    other.mapping_ = nullptr;
    other.mapping_bytes_ = 0;
    other.records_ = nullptr;
    other.count_ = 0;
  }
  return *this;
}

FactorBase FactorBase::Build(factor_t p_max) {
  FactorBase fb;
  fb.p_max_ = p_max;
  fb.owned_.reserve(EstimateSplitPrimes(p_max));

  const SegmentedSieve sieve(p_max);
  sieve.for_each_split_prime([&](factor_t p) {
    const Gauss z = decompose_prime(p); // throws if p is not split
    fb.owned_.push_back({static_cast<std::uint32_t>(z.re),
                         static_cast<std::uint32_t>(z.im)});
  });
  fb.records_ = fb.owned_.data();
  fb.count_ = fb.owned_.size();
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
    buffer.push_back(
        {static_cast<std::uint32_t>(z.re), static_cast<std::uint32_t>(z.im)});
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
  std::rewind(in.get());

  const std::size_t count =
      RecordCount(static_cast<std::size_t>(bytes), path);

  FactorBase fb;
  // Read in blocks and stop at the bound, rather than slurping the file and
  // truncating afterwards. The difference is not subtle: the shipped artifact
  // is 2.75 GB, of which a B = 10^10 campaign needs the first 34 MB, and the
  // naive version allocated all of it.
  std::vector<Record> block(1u << 16);
  std::size_t remaining = count;
  factor_t previous = 0;
  bool done = false;

  while (remaining > 0 && !done) {
    const std::size_t want = std::min(remaining, block.size());
    if (std::fread(block.data(), sizeof(Record), want, in.get()) != want)
      throw std::runtime_error("factorbase: short read from " + path);
    remaining -= want;

    for (std::size_t i = 0; i < want; ++i) {
      const Record &r = block[i];
      if (r.a < r.b || r.b == 0) {
        throw std::runtime_error("factorbase: " + path + " record " +
                                 std::to_string(fb.owned_.size()) +
                                 " violates a >= b > 0");
      }
      const factor_t p = PrimeOf(r);
      if (p <= previous) {
        throw std::runtime_error("factorbase: " + path +
                                 " is not ascending at p = " +
                                 std::to_string(p));
      }
      if (p_max != 0 && p > p_max) {
        done = true;
        break;
      }
      previous = p;
      fb.owned_.push_back(r);
    }
  }
  if (p_max != 0 && !done && previous < p_max) {
    throw std::runtime_error("factorbase: " + path + " reaches only " +
                             std::to_string(previous) +
                             ", short of the requested bound " +
                             std::to_string(p_max));
  }
  const std::size_t used = fb.owned_.size();
  fb.owned_.shrink_to_fit();
  fb.records_ = fb.owned_.data();
  fb.count_ = used;
  fb.p_max_ = p_max != 0 ? p_max
                         : (used == 0 ? 0 : PrimeOf(fb.owned_[used - 1]));
  return fb;
}

FactorBase FactorBase::Map(const std::string &path, factor_t p_max) {
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw std::runtime_error("factorbase: cannot open " + path);

  struct stat st {};
  if (::fstat(fd, &st) != 0) {
    ::close(fd);
    throw std::runtime_error("factorbase: cannot stat " + path);
  }
  const auto bytes = static_cast<std::size_t>(st.st_size);
  std::size_t count = 0;
  try {
    count = RecordCount(bytes, path);
  } catch (...) {
    ::close(fd);
    throw;
  }

  FactorBase fb;
  if (count == 0) {
    ::close(fd);
    fb.p_max_ = p_max;
    return fb;
  }

  void *addr = ::mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
  // The descriptor is not needed once the mapping exists; the mapping keeps the
  // underlying file alive by itself.
  ::close(fd);
  if (addr == MAP_FAILED)
    throw std::runtime_error("factorbase: cannot map " + path);

  fb.mapping_ = addr;
  fb.mapping_bytes_ = bytes;
  const auto *records = static_cast<const Record *>(addr);
  fb.records_ = records;
  fb.count_ = count;

  // **Validation here is sampled, not exhaustive, and that is the entire point
  // of mapping.** Walking all 343M records to check ordering would fault in the
  // whole 2.7 GB file, at which point Map costs exactly what Load costs and
  // buys only evictability. Measured: the exhaustive version pinned 2,623 MB
  // resident, identical to reading it.
  //
  // What the sample catches is what actually goes wrong with an artifact:
  // truncation, misalignment, a wrong or partially written file. What it cannot
  // catch is a single flipped record in the middle -- and that is covered far
  // better by cpu/audit/, which re-derives decompositions with decompose_prime
  // rather than merely checking that the file is sorted. validate_all() is
  // available for a caller that wants the full pass anyway.
  try {
    fb.ValidateSample(path);
    // Binary search for the cut, O(log n) page touches rather than O(n).
    fb.count_ = p_max == 0 ? count : fb.upper_bound(p_max);
    if (fb.count_ == 0) {
      throw std::runtime_error("factorbase: " + path +
                               " has no primes at or below " +
                               std::to_string(p_max));
    }
    if (p_max != 0 && fb.count_ == count &&
        PrimeOf(records[count - 1]) < p_max) {
      throw std::runtime_error(
          "factorbase: " + path + " reaches only " +
          std::to_string(PrimeOf(records[count - 1])) +
          ", short of the requested bound " + std::to_string(p_max));
    }
  } catch (...) {
    fb.Reset();
    throw;
  }

  fb.p_max_ = p_max != 0 ? p_max : PrimeOf(records[fb.count_ - 1]);
  return fb;
}

void FactorBase::ValidateSample(const std::string &path) const {
  if (count_ == 0)
    return;
  // Head and tail in full -- a truncated or misaligned file shows up there --
  // plus evenly spaced probes checking that the ordering holds globally.
  const std::size_t edge = std::min<std::size_t>(count_, 1024);
  const auto check_run = [&](std::size_t from, std::size_t to) {
    factor_t previous = from == 0 ? 0 : PrimeOf(records_[from - 1]);
    for (std::size_t i = from; i < to; ++i) {
      const Record &r = records_[i];
      if (r.a < r.b || r.b == 0) {
        throw std::runtime_error("factorbase: " + path + " record " +
                                 std::to_string(i) + " violates a >= b > 0");
      }
      const factor_t p = PrimeOf(r);
      if (p <= previous) {
        throw std::runtime_error("factorbase: " + path +
                                 " is not ascending at record " +
                                 std::to_string(i));
      }
      previous = p;
    }
  };
  check_run(0, edge);
  check_run(count_ - std::min(count_, edge), count_);

  constexpr std::size_t kProbes = 4096;
  const std::size_t stride = std::max<std::size_t>(1, count_ / kProbes);
  factor_t previous = 0;
  for (std::size_t i = 0; i < count_; i += stride) {
    const Record &r = records_[i];
    if (r.a < r.b || r.b == 0) {
      throw std::runtime_error("factorbase: " + path + " record " +
                               std::to_string(i) + " violates a >= b > 0");
    }
    const factor_t p = PrimeOf(r);
    if (p <= previous) {
      throw std::runtime_error("factorbase: " + path +
                               " is not ascending at record " +
                               std::to_string(i));
    }
    previous = p;
  }
}

void FactorBase::validate_all(const std::string &path) const {
  ValidatePrefix(records_, count_, 0, path);
}

std::vector<factor_t> FactorBase::primes_vector() const {
  std::vector<factor_t> out;
  out.reserve(count_);
  for (std::size_t i = 0; i < count_; ++i)
    out.push_back(prime(static_cast<factor_index_t>(i)));
  return out;
}

std::size_t FactorBase::upper_bound(factor_t limit) const {
  std::size_t lo = 0;
  std::size_t hi = count_;
  while (lo < hi) {
    const std::size_t mid = lo + (hi - lo) / 2;
    if (prime(static_cast<factor_index_t>(mid)) <= limit)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

} // namespace pte::cpu
