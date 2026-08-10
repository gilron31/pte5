// Fixed-width unsigned integers of any limb count, shared by the host and the
// device. One thread owns one value; there is no cooperation between lanes.
//
// ---------------------------------------------------------------------------
// What this is for
// ---------------------------------------------------------------------------
//
// The exhaustive campaign never needs it: at B = 10^12 a coordinate is 20 bits,
// nu is 40, and the exact pair key is 81 -- all of which cpu/core/precision.h
// covers with int32/int64/__int128. **S-E, the champion meet, is outside that
// policy**, and this file is the width it needs
// (drafts/gpu_implementation_plan.md 3.3):
//
//   quantity                    champion (L1 ~ 2.5e36)   type here
//   ---------------------------------------------------------------
//   Gauss coordinates m, n      61 bits                  uint32 x 2
//   nu = m^2 - n^2              122 bits                 U128 = UInt<4>
//   pair key nu_a^2 + nu_b^2    ~244 bits                U256 = UInt<8>
//
// Two properties are load-bearing and neither is negotiable:
//
//   1. **The windowing comparison must be exact.** S-E partitions the pair
//      space by key *value*; a comparison that is off by one ulp puts a pair in
//      the wrong window, where it is silently never compared against its
//      partner. That is a lost solution with no diagnostic -- the one failure
//      mode Part IV of the design doc says actually matters. Hence exact
//      integer compares here, never `double`.
//   2. **The layout is the host/device contract.** A value is an array of
//      32-bit limbs, least significant first, and it crosses the PCIe bus as
//      raw bytes. It is a trivially copyable aggregate with no padding, so
//      sizeof(UInt<N>) == 4N on both sides; the static_asserts at the bottom
//      say so out loud.
//
// ---------------------------------------------------------------------------
// Why 32-bit limbs, and why the carry chain is inline PTX
// ---------------------------------------------------------------------------
//
// 32-bit limbs are the convention cpu/bigint/bigint.h and gpu/bigint/bigint.h
// already state ("a fixed-width array of 32-bit limbs, least-significant limb
// first"), so `pte::cpu::to_limbs` / `from_limbs` bridge to GMP with no
// repacking, and gpu/bigint's schoolbook kernel speaks the same layout.
//
// They are also what the hardware is. **sm_89 has no 64-bit integer ALU**: a
// 64-bit add is emitted as IADD3 + IADD3.X, two 32-bit adds through the carry
// flag. 64-bit limbs are not a shortcut past that, only a less honest
// description of it -- and the description turns out to cost real time, because
// what the compiler does with the *carry* depends on how the limbs are spelled.
//
// Measured, and the SASS read (bench/gpu P8; SASS counted on an add-only kernel
// so the histogram is the add and nothing else):
//
//   256-bit add           SASS ops   op-pair/s   relative   S-E's 1.374e11
//   ------------------------------------------------------------------------
//   4 x uint64, C           ~25       1.14e11      1.00         1.21 s
//   8 x uint32, C            14       1.70e11      1.49         0.81 s
//   8 x uint32, inline PTX    7       2.53e11      2.22         0.54 s
//
// (All three rows at 2.12 GHz in that run, so the ratios are direct. Read the
// clock the harness prints beside any rerun: a shorter kernel boosts higher,
// and an earlier run at 1.89/2.25/2.47 GHz gave the same 2.22x once
// normalised.)
//
// The instruction counts explain the times almost exactly, and they are the
// real finding: **with 64-bit limbs the compiler stops using the carry flag.**
// It emits the add, then rebuilds the carry with ISETP + SEL per limb -- 3 SEL
// and 10 ISETP for one 256-bit add. With 32-bit limbs and a uint64 accumulate
// it recovers the chain on its own (8 IADD3 + 6 IADD3.X); with the explicit
// `add.cc.u32` / `addc.cc.u32` / `addc.u32` sequence it is the chain and
// nothing else, one instruction per limb.
//
// So the PTX is in for the instruction count, not for the half second: against
// P4's 138 s of collision detection, S-E's arithmetic is under 1% either way.
// It matters wherever the arithmetic *is* the loop, and it costs nothing to
// have been written correctly once.
//
// Deliberately minimal: compare, add, subtract, shift, and full-width multiply.
// That is exactly what S-E uses. There is no division and no signed form.
#pragma once

#include "core/hostdev.h"

#include <cstddef>
#include <cstdint>

namespace pte::core {

// ---------------------------------------------------------------------------
// The type
// ---------------------------------------------------------------------------

// An unsigned integer of exactly kLimbs 32-bit limbs, least significant first.
//
// An aggregate with a public array and no constructors, on purpose: it has to
// be trivially copyable to be a kernel argument and to cross the bus by
// memcpy, and it has to be brace-initialisable to be a compile-time constant.
//
// **Over-aligned to 16 bytes when the width allows it.** The layout is
// unchanged -- still kLimbs consecutive 32-bit limbs, still 4*kLimbs bytes,
// still memcpy-able to the device -- but alignment is what lets nvcc vectorize
// a load of one of these. Confirmed in the shipped kernel's SASS:
//
//   alignof 4    8 x LDG.E      (32-bit) per U256, each a stride-32B gather
//   alignof 16   2 x LDG.E.128  per U256
//
// **And it makes no measurable difference**, which is the more useful half of
// the result. Sweeping the champion kernel's in-kernel rate across four table
// sizes, before and after: 2.76 / 2.22 / 1.05 / 0.964 e9 pair/s against
// 2.74 / 2.19 / 1.05 / 0.962. Identical inside the noise. The sq[] loads are
// not what the kernel waits on -- the hash-table atomics are (see
// gpu/champion/champion_impl.cu), and the same sweep shows *those* moving 2.9x
// with the table's address footprint. Kept anyway because it is strictly fewer
// instructions for no cost, and because it will matter if the loads ever stop
// being free; not kept under any claim that it sped anything up.
//
// A first version of this file asserted `alignof(U256) == 4` on the theory that
// minimal alignment protected the raw-bytes contract. That theory was wrong:
// the contract is size and limb order, and alignment touches neither. The
// condition keeps it honest at every width -- 16-byte alignment is applied only
// when 4*kLimbs is already a multiple of 16, so sizeof never grows (a UInt<6>
// would otherwise pad 24 bytes to 32).
template <int kLimbs>
struct alignas(kLimbs % 4 == 0 ? 16 : alignof(std::uint32_t)) UInt {
  static_assert(kLimbs >= 2, "UInt needs at least two limbs: the carry chain "
                             "is written as first / middle / last");
  static constexpr int kBits = 32 * kLimbs;

  std::uint32_t w[kLimbs]{};
};

// The two widths S-E actually instantiates.
using U128 = UInt<4>; // a champion nu (122 bits)
using U256 = UInt<8>; // a champion pair key (~244 bits)

template <int N> PTE_HOST_DEVICE inline UInt<N> uint_from(std::uint64_t v) {
  UInt<N> r;
  r.w[0] = static_cast<std::uint32_t>(v);
  r.w[1] = static_cast<std::uint32_t>(v >> 32);
  for (int i = 2; i < N; ++i)
    r.w[i] = 0;
  return r;
}

// Widen, zero-extending. Compile-time checked so a narrowing typo cannot
// silently truncate a key.
template <int M, int N> PTE_HOST_DEVICE inline UInt<M> widen(const UInt<N> &a) {
  static_assert(M >= N, "widen must not truncate");
  UInt<M> r;
  for (int i = 0; i < N; ++i)
    r.w[i] = a.w[i];
  for (int i = N; i < M; ++i)
    r.w[i] = 0;
  return r;
}

PTE_HOST_DEVICE inline U128 u128_from(std::uint64_t lo, std::uint64_t hi = 0) {
  U128 r;
  r.w[0] = static_cast<std::uint32_t>(lo);
  r.w[1] = static_cast<std::uint32_t>(lo >> 32);
  r.w[2] = static_cast<std::uint32_t>(hi);
  r.w[3] = static_cast<std::uint32_t>(hi >> 32);
  return r;
}

PTE_HOST_DEVICE inline U256 u256_from(std::uint64_t v) {
  return uint_from<8>(v);
}
PTE_HOST_DEVICE inline U256 u256_from(const U128 &a) { return widen<8, 4>(a); }

// The i'th 64-bit chunk, for callers that want to hash or print a value
// without knowing the limb width. Little-endian, like the limbs.
template <int N>
PTE_HOST_DEVICE inline std::uint64_t chunk64(const UInt<N> &a, int i) {
  return static_cast<std::uint64_t>(a.w[2 * i]) |
         (static_cast<std::uint64_t>(a.w[2 * i + 1]) << 32);
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------
//
// Returns -1, 0 or +1. Lexicographic from the most significant limb down,
// which for a little-endian limb array means counting backwards -- the loop
// direction is the whole content of the function, and getting it wrong compares
// two numbers by their low bits, which is exactly the silent misordering
// property 1 above is about.

template <int N>
PTE_HOST_DEVICE inline int cmp(const UInt<N> &a, const UInt<N> &b) {
  for (int i = N - 1; i >= 0; --i) {
    if (a.w[i] != b.w[i])
      return a.w[i] < b.w[i] ? -1 : 1;
  }
  return 0;
}

template <int N>
PTE_HOST_DEVICE inline bool operator==(const UInt<N> &a, const UInt<N> &b) {
  return cmp(a, b) == 0;
}
template <int N>
PTE_HOST_DEVICE inline bool operator!=(const UInt<N> &a, const UInt<N> &b) {
  return cmp(a, b) != 0;
}
template <int N>
PTE_HOST_DEVICE inline bool operator<(const UInt<N> &a, const UInt<N> &b) {
  return cmp(a, b) < 0;
}

template <int N> PTE_HOST_DEVICE inline bool is_zero(const UInt<N> &a) {
  std::uint32_t acc = 0;
  for (int i = 0; i < N; ++i)
    acc |= a.w[i];
  return acc == 0;
}

// ---------------------------------------------------------------------------
// Addition and subtraction
// ---------------------------------------------------------------------------
//
// Both wrap modulo 2^(32N) and neither reports the carry out. That is safe for
// every use in S-E and it is worth saying why rather than trusting it: a
// champion key is nu_a^2 + nu_b^2 <= 2 * L1^2, and L1 ~ 2.5e36 < 2^122 gives
// 2^245 -- eleven bits of headroom in a U256. `sub` is only ever called on
// lo <= hi window bounds. A caller outside that domain gets modular arithmetic,
// which is why neither is spelled `operator+` on a type that also compares.
//
// Two implementations, and they must agree bit for bit (core/hostdev.h; checked
// directly in gpu/bigint/uint_device_test.cu). See the table at the top of the
// file for what the difference buys and what it does not.
//
//   host    accumulate in uint64 and shift the carry out. Portable, and x86
//           has ADC anyway; the host path is a reference, not a hot loop.
//   device  the explicit carry chain, one instruction per limb.

template <int N>
PTE_HOST_DEVICE inline UInt<N> add(const UInt<N> &a, const UInt<N> &b) {
  UInt<N> r;
#if defined(__CUDA_ARCH__)
  // asm volatile, not plain asm: the carry flag is an implicit dependency the
  // compiler cannot see, so nothing but the volatile ordering keeps these
  // instructions adjacent. This is the CGBN idiom.
  asm volatile("add.cc.u32 %0, %1, %2;"
               : "=r"(r.w[0])
               : "r"(a.w[0]), "r"(b.w[0]));
#pragma unroll
  for (int i = 1; i < N - 1; ++i)
    asm volatile("addc.cc.u32 %0, %1, %2;"
                 : "=r"(r.w[i])
                 : "r"(a.w[i]), "r"(b.w[i]));
  asm volatile("addc.u32 %0, %1, %2;"
               : "=r"(r.w[N - 1])
               : "r"(a.w[N - 1]), "r"(b.w[N - 1]));
#else
  std::uint64_t carry = 0;
  for (int i = 0; i < N; ++i) {
    const std::uint64_t s = static_cast<std::uint64_t>(a.w[i]) + b.w[i] + carry;
    r.w[i] = static_cast<std::uint32_t>(s);
    carry = s >> 32;
  }
#endif
  return r;
}

template <int N>
PTE_HOST_DEVICE inline UInt<N> sub(const UInt<N> &a, const UInt<N> &b) {
  UInt<N> r;
#if defined(__CUDA_ARCH__)
  asm volatile("sub.cc.u32 %0, %1, %2;"
               : "=r"(r.w[0])
               : "r"(a.w[0]), "r"(b.w[0]));
#pragma unroll
  for (int i = 1; i < N - 1; ++i)
    asm volatile("subc.cc.u32 %0, %1, %2;"
                 : "=r"(r.w[i])
                 : "r"(a.w[i]), "r"(b.w[i]));
  asm volatile("subc.u32 %0, %1, %2;"
               : "=r"(r.w[N - 1])
               : "r"(a.w[N - 1]), "r"(b.w[N - 1]));
#else
  std::uint64_t borrow = 0;
  for (int i = 0; i < N; ++i) {
    // Wraps by design when a < b; the borrow is then the high half being all
    // ones, which is what the shift reads off.
    const std::uint64_t d =
        static_cast<std::uint64_t>(a.w[i]) - b.w[i] - borrow;
    r.w[i] = static_cast<std::uint32_t>(d);
    borrow = (d >> 32) & 1u;
  }
#endif
  return r;
}

// ---------------------------------------------------------------------------
// Shifts
// ---------------------------------------------------------------------------
//
// Defined for 0 <= n < kBits; anything else returns zero rather than invoking
// the undefined behaviour of a 32-bit shift by 32. The explicit `bits == 0`
// arm is not redundant: `x >> 32` is UB, not zero, and that is the one input a
// shift helper is usually wrong for.

template <int N> PTE_HOST_DEVICE inline UInt<N> shl(const UInt<N> &a, int n) {
  UInt<N> r = UInt<N>{};
  if (n < 0 || n >= 32 * N)
    return r;
  const int words = n / 32, bits = n % 32;
  for (int i = N - 1; i >= words; --i) {
    const int src = i - words;
    std::uint32_t v = a.w[src] << bits;
    if (bits != 0 && src > 0)
      v |= a.w[src - 1] >> (32 - bits);
    r.w[i] = v;
  }
  return r;
}

template <int N> PTE_HOST_DEVICE inline UInt<N> shr(const UInt<N> &a, int n) {
  UInt<N> r = UInt<N>{};
  if (n < 0 || n >= 32 * N)
    return r;
  const int words = n / 32, bits = n % 32;
  for (int i = 0; i + words < N; ++i) {
    const int src = i + words;
    std::uint32_t v = a.w[src] >> bits;
    if (bits != 0 && src + 1 < N)
      v |= a.w[src + 1] << (32 - bits);
    r.w[i] = v;
  }
  return r;
}

// Position of the highest set bit plus one; 0 for zero.
template <int N> PTE_HOST_DEVICE inline int bit_length(const UInt<N> &a) {
  for (int i = N - 1; i >= 0; --i) {
    if (a.w[i] == 0)
      continue;
    int b = 0;
    for (std::uint32_t v = a.w[i]; v != 0; v >>= 1)
      ++b;
    return i * 32 + b;
  }
  return 0;
}

// ---------------------------------------------------------------------------
// Multiplication
// ---------------------------------------------------------------------------
//
// Full width: N limbs times N limbs into 2N, so nothing is discarded and there
// is no overflow to reason about. Schoolbook, one thread, N^2 limb products.
//
// Left as portable C rather than PTX like the carry chain above, and the reason
// is the workload rather than taste: S-E multiplies exactly D ~ 2^19 times (one
// square per point, once per run) against 1.374e11 additions in the pair loop.
// Optimizing it would be optimizing 0.0000004% of the search.
//
// The accumulator cannot overflow, which is the one thing to check in a
// schoolbook loop: (2^32-1)^2 + 2*(2^32-1) = 2^64 - 1 exactly.
template <int N>
PTE_HOST_DEVICE inline UInt<2 * N> mul(const UInt<N> &a, const UInt<N> &b) {
  UInt<2 * N> r = UInt<2 * N>{};
  for (int i = 0; i < N; ++i) {
    std::uint32_t carry = 0;
    for (int j = 0; j < N; ++j) {
      const std::uint64_t t =
          static_cast<std::uint64_t>(a.w[i]) * b.w[j] + r.w[i + j] + carry;
      r.w[i + j] = static_cast<std::uint32_t>(t);
      carry = static_cast<std::uint32_t>(t >> 32);
    }
    // Safe as an assignment rather than an add: row i writes limbs i..i+N-1,
    // so limb i+N has not been touched by any earlier row.
    r.w[i + N] = carry;
  }
  return r;
}

// nu^2, the only multiplication S-E performs. Named because `mul(x, x)` at a
// call site reads like a typo.
template <int N> PTE_HOST_DEVICE inline UInt<2 * N> square(const UInt<N> &a) {
  return mul(a, a);
}

// ---------------------------------------------------------------------------
// The layout contract
// ---------------------------------------------------------------------------

// Size and limb order are the contract; alignment is free to be whatever makes
// the device fast, as long as it does not disturb them.
static_assert(sizeof(U128) == 16, "U128 must be exactly four 32-bit limbs");
static_assert(sizeof(U256) == 32, "U256 must be exactly eight 32-bit limbs");
static_assert(sizeof(UInt<6>) == 24,
              "over-alignment must never pad a width it does not divide");
static_assert(alignof(U256) == 16,
              "16-byte alignment is what lets a U256 load vectorize to "
              "LDG.E.128; see the note on the declaration");
static_assert(offsetof(U256, w) == 0,
              "the limb array must start at offset 0: the device reads these "
              "bytes as an array of uint32");

} // namespace pte::core
