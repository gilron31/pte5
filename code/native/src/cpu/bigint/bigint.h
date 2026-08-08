// Big-integer arithmetic on the CPU, backed by GMP.
//
// Also provides the limb export/import used to hand values to the GPU kernels:
// a fixed-width array of 32-bit limbs, least-significant limb first. That
// layout is the contract shared with src/gpu/gpu.h.
#ifndef PTE_CPU_BIGINT_H_
#define PTE_CPU_BIGINT_H_

#include <gmpxx.h>

#include <cstddef>
#include <cstdint>
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

} // namespace pte::cpu

#endif // PTE_CPU_BIGINT_H_
