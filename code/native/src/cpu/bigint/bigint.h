// Big-integer arithmetic on the CPU, backed by GMP.
//
// Also provides the limb export/import used to hand values to the GPU kernels:
// a fixed-width array of 32-bit limbs, least-significant limb first. That
// layout is the contract shared with src/gpu/bigint/bigint.h.
#ifndef PTE_CPU_BIGINT_H_
#define PTE_CPU_BIGINT_H_

#include "core/uint.h"

#include <gmpxx.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace pte::cpu {

// Number of bits carried by one limb of the GPU-facing representation.
inline constexpr int kLimbBits = 32;

mpz_class add(const mpz_class &a, const mpz_class &b);
mpz_class mul(const mpz_class &a, const mpz_class &b);

// Decimal-string convenience wrappers, handy from a REPL or a demo.
std::string add_dec(const std::string &a, const std::string &b);
std::string mul_dec(const std::string &a, const std::string &b);

// Number of 32-bit limbs needed to hold |z|. Zero needs zero limbs.
std::size_t limb_count(const mpz_class &z);

// Writes |z| into out[0..nlimbs), least-significant limb first, zero-padded.
// Throws std::invalid_argument if z is negative or does not fit in nlimbs.
void to_limbs(const mpz_class &z, std::uint32_t *out, std::size_t nlimbs);

// Convenience overload returning a freshly sized vector.
std::vector<std::uint32_t> to_limbs(const mpz_class &z, std::size_t nlimbs);

// Inverse of to_limbs.
mpz_class from_limbs(const std::uint32_t *in, std::size_t nlimbs);

// ---------------------------------------------------------------------------
// The bridge to the fixed-width type
// ---------------------------------------------------------------------------
//
// core::UInt<N> *is* the limb array above -- N 32-bit limbs, least significant
// first -- so these two are a cast with a bounds check rather than a
// conversion. That is the payoff of keeping the shared type on the documented
// convention instead of on 64-bit limbs: no repacking anywhere, and GMP stays
// on the host side of the seam where cpu/core/precision.h wants it (layer 3,
// the certificate layer).
//
// S-E uses both directions: champion point sets are built in mpz_class because
// the coordinates are 61 bits and Cornacchia is a GMP job, then narrowed to
// core::U128 for the meet; a confirmed hit goes back to mpz_class to be
// certified.
template <int N> core::UInt<N> to_uint(const mpz_class &z) {
  core::UInt<N> r{};
  to_limbs(z, r.w, static_cast<std::size_t>(N)); // throws if it does not fit
  return r;
}

template <int N> mpz_class from_uint(const core::UInt<N> &v) {
  return from_limbs(v.w, static_cast<std::size_t>(N));
}

} // namespace pte::cpu

#endif // PTE_CPU_BIGINT_H_
