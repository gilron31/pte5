#include "bigint.h"

#include <stdexcept>

namespace pte::cpu {
namespace {

// mpz_export/mpz_import word format shared by to_limbs and from_limbs:
//   order = -1  -> least significant word first
//   size  =  4  -> 32-bit words
//   endian = 0  -> native byte order within each word
//   nails  = 0  -> use all bits
constexpr int kOrder = -1;
constexpr std::size_t kWordBytes = 4;
constexpr int kEndian = 0;
constexpr std::size_t kNails = 0;

}  // namespace

mpz_class add(const mpz_class& a, const mpz_class& b) { return a + b; }

mpz_class mul(const mpz_class& a, const mpz_class& b) { return a * b; }

std::string add_dec(const std::string& a, const std::string& b) {
  return add(mpz_class(a), mpz_class(b)).get_str();
}

std::string mul_dec(const std::string& a, const std::string& b) {
  return mul(mpz_class(a), mpz_class(b)).get_str();
}

std::size_t limb_count(const mpz_class& z) {
  if (z == 0) return 0;
  // Ceiling division of the bit length by the limb width.
  const std::size_t bits = mpz_sizeinbase(z.get_mpz_t(), 2);
  return (bits + kLimbBits - 1) / kLimbBits;
}

void to_limbs(const mpz_class& z, std::uint32_t* out, std::size_t nlimbs) {
  if (z < 0) {
    throw std::invalid_argument("to_limbs: negative values are not representable");
  }
  const std::size_t needed = limb_count(z);
  if (needed > nlimbs) {
    throw std::invalid_argument("to_limbs: value needs " + std::to_string(needed) +
                                " limbs but only " + std::to_string(nlimbs) +
                                " were provided");
  }
  // mpz_export writes nothing at all when z == 0, so the zero fill must happen
  // after the export rather than relying on it to touch every limb.
  std::size_t written = 0;
  if (needed > 0) {
    mpz_export(out, &written, kOrder, kWordBytes, kEndian, kNails, z.get_mpz_t());
  }
  for (std::size_t i = written; i < nlimbs; ++i) out[i] = 0;
}

std::vector<std::uint32_t> to_limbs(const mpz_class& z, std::size_t nlimbs) {
  std::vector<std::uint32_t> limbs(nlimbs);
  to_limbs(z, limbs.data(), nlimbs);
  return limbs;
}

mpz_class from_limbs(const std::uint32_t* in, std::size_t nlimbs) {
  mpz_class z;
  mpz_import(z.get_mpz_t(), nlimbs, kOrder, kWordBytes, kEndian, kNails, in);
  return z;
}

}  // namespace pte::cpu
