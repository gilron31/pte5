# Code review — responses and TODO

Working file for the review of `code/native`. Part A answers the inline
`GILRO:` comments one by one; Part B is my own audit of the same code, which
turned up things the inline comments did not; Part C is the consolidated,
prioritized TODO.

Severity labels: **[S]** silent-wrong-answer risk (the failure mode that
actually threatens an exclusion claim), **[C]** correctness/robustness,
**[P]** performance, **[A]** architecture/hygiene.

No code has been changed. Every claim below that has a number attached was
checked; the checks are noted inline.

---

## Part A — responses to the inline comments

### A1. `DeviceBuffer` belongs in `device_buffer.h`, not `check.cuh` — **agreed**

`src/gpu/check.cuh:34`. Correct, and it is the cheapest item on the list. Error
checking and an owning allocation are unrelated concerns that got stapled
together because the scaffold only had two kernels.

Split into `check.cuh` (the macro + `ThrowIfCudaError`) and
`device_buffer.cuh` (the RAII type). Do it in the same pass as A4 so the
directory move happens once.

### A2. `HostBuffer<T>` with pinned memory — **agreed on direction, two changes**

`src/gpu/check.cuh:72`. The reasoning is right: pinned memory is what unlocks
`cudaMemcpyAsync`, and pageable→pinned staging costs an extra copy that the
driver is doing invisibly today anyway.

Two things I would change about the proposal:

**(a) The producer must write *directly into* the `HostBuffer`, not copy into
it.** The comment describes "each cpu input should be first moved into a
preallocated HostBuffer, and only then to a DeviceBuffer". If that move is a
`std::vector` → `HostBuffer` copy, we have *added* a full pass over the data to
save a driver-internal one. For the enumerator the fix is easy — the CSR packing
loop in `enum_demo.cc` should emit into the `HostBuffer` in the first place — but
the API has to make that the natural thing to do, e.g. `HostBuffer` exposing a
`std::span` the producer fills.

**(b) Pinned allocation is slow and a limited system resource** —
`cudaHostAlloc` is far more expensive than `malloc` and pins physical pages. It
must be allocated once and reused across batches, never per call. That pushes
toward a session/arena object owning the buffers, which is the same conclusion
A4 reaches from a different direction.

**Strongly agree on "build a test or example that demonstrates why we need
HostBuffer".** That is the right way to settle it, and it is half a day: a
benchmark that copies the same buffer pageable vs. pinned, sync vs. async, and
reports achieved GB/s. Worth doing *before* the refactor, because it also tells
us whether the transfer path is worth optimizing at all given that the plan is
to stop transferring norms entirely (see P2 below).

### A3. `#pragma once` instead of include guards — **agreed, with one caveat**

`src/gpu/gpu.h:11`. Fine by me. Caveat: all five existing headers use guards, so
this should be one mechanical sweep across the whole tree rather than a new
convention that half the files follow. `#pragma once` is not ISO standard but is
supported by GCC, Clang, MSVC and NVCC, which is the entire universe this
project targets.

### A4/A5. Break up the monolithic `gpu.h` into `app.h` / `app_impl.h` / `app_impl.cu` — **agreed on the split, pushing back on the middle layer's signature**

`src/gpu/gpu.h:16` and `:38`. The problem is real and getting worse: `gpu.h` now
carries three unrelated APIs (bigint, matmul, meet), and `meet_batch`'s
declaration sits between the matmul and the limb-layout comment with nothing
connecting them. Per-app files are right, and tests going through the `app.h`
path is right.

Where I want to push back, concretely:

**(a) "app.h receives std::vector, allocates HostBuffer and DeviceBuffer,
performs the copies, and invokes app_impl" bakes per-call allocation into the
API.** That is exactly what `meet_batch` does today, and it is why the demo
spends its time in `cudaMalloc`/`cudaMemcpy` rather than in the kernel. For the
campaign we need buffers that persist across thousands of batches. If the public
entry point owns allocation, every caller pays it every time and there is no
seam at which to hoist it.

Suggested shape instead — same three layers, but ownership hoisted out:

```
meet.h        plain C++. Two entry points:
                - a convenience one taking spans (tests, demo)
                - one taking a MeetSession& that owns the buffers (campaign)
meet_session  opaque handle; definition is CUDA-side (pimpl), so meet.h stays
              free of <cuda_runtime.h>
meet_impl.cuh device-side signatures over device pointers
meet_impl.cu  the kernel
```

**(b) There is a seam conflict in the proposal that needs deciding explicitly.**
`app.h` cannot name `HostBuffer` or `DeviceBuffer` in a signature if those are
templates defined in a `.cuh` that includes `<cuda_runtime.h>` — the demo and the
pybind module are compiled by `g++` and including that header breaks the
build. Options: (i) pimpl/opaque handle as above; (ii) `app.h` becomes
CUDA-aware and a thinner `app_pub.h` keeps the plain-C++ promise; (iii) drop the
seam rule. I would take (i) — it keeps the property that made the GPU tests
provably independent of `pte_cpu`.

### A6. Named types, and `<cstdint>` instead of `unsigned long long` — **agreed, and this is hiding a real bug**

`src/gpu/meet_kernel.cu:27`. Agreed on its own terms, but it is more than style:
the type sprawl is currently papered over with `reinterpret_cast`, and one of
those casts hides a genuine scaling bug (see **C1** below).

Today the same quantity is spelled three ways — `std::uint64_t` in `gpu.h`,
`unsigned long long` in the kernel, bridged by
`reinterpret_cast<const unsigned long long *>(nu)`; and offsets are
`std::int32_t` publicly but `int` in the kernel. Those casts are also
formally type-punning between distinct types that merely happen to share a
representation.

Proposed, in a header shared by host and device:

```cpp
using nu_t    = std::uint64_t;   // point invariant; wraps deliberately
using key_t   = std::uint64_t;   // pair key, mod 2^64
using index_t = std::int64_t;    // CSR offsets  <- widened, see C1
using count_t = std::uint64_t;   // hit counter  <- widened, see C1
```

One constraint to respect: `atomicAdd` requires exactly
`unsigned long long int` for the 64-bit overload, so `count_t` needs a
device-side alias that is spelled that way even if the host-side name differs.

### A7. Make the output buffer a class with its own tests — **agreed on substance, one correction**

`src/gpu/meet_kernel.cu:39`. Agreed, and "this C-type programming is very error
prone" is fair — `Emit` takes eight parameters, three of which (`hits`,
`max_hits`, `hit_count`) are one object smeared across three arguments, and the
overflow contract lives half in the kernel and half in `meet_batch`.

The correction: **it is not a hash table.** It is an append-only output buffer
with an atomic bump allocator. Naming it a hash table would mislead the next
reader, and it also matters for planning, because a *real* hash table is a
separate thing we will need later — the algorithm doc's A4 switches from
direct comparison to hashing for `D >= 10`. Two different components:

- now: `DeviceAppendBuffer<T>` — `(ptr, capacity, atomic counter)`, a
  `try_push` returning whether it fit, and the "count keeps rising past
  capacity so the caller knows the retry size" semantics as an explicit,
  documented, independently tested contract. The truncation behaviour is
  currently only tested indirectly through `ReportsTruncationWithoutOverrunning`.
- later: an actual device hash table for the large-`D` bucket.

### A8. Amortize the squarings — **agreed, and it is worth more than it looks**

`src/gpu/meet_kernel.cu:77`. Agreed, and I would not count on the compiler here:
`nu` and `hits` are both pointers into global memory with no `__restrict__`, so
nvcc cannot prove `Emit`'s writes do not alias `v[]`, and it is therefore not
free to hoist the loads or the squares out of the loop.

Counting the innermost body as written — `kij` and `kik` are already hoisted, so
per 4-subset we do `PairKey` four times = **8 multiplies + 4 adds**. Precomputing
`sq[t] = v[t] * v[t]` once per norm turns every key into a single add: **0
multiplies + 4 adds**, for `D` multiplies of setup. That is a 3–4× cut in the
hot loop's ALU work.

Two caveats worth writing down before implementing:

- Storage. `sq[]` is `D` × 8 bytes per *thread*. At `D <= 8` (97% of norms) it
  is 64 bytes and lives in registers; at the `B = 10^12` maximum of `D = 192` it
  is 1.5 KB per thread, which will spill. So this lands naturally together with
  the `D`-bucketing work, not before it.
- `const __restrict__` on the `nu` parameter is a one-line change that should go
  in regardless, and may recover part of this on its own.

---

## Part B — my audit of the same code

Things I found reviewing my own code that the inline comments did not cover.
The first three are the ones I would fix before anything in Part A.

### C1 [S]. `index_t` is 32-bit — the campaign overflows it silently

`meet_batch` takes `const std::int32_t *offsets` and returns `int`; the kernel
uses `int` internally and the hit counter is `unsigned int`.

- Total points in a batch is capped at 2^31. The extrapolation in the algorithm
  doc is ~10^11 points at `B = 10^12`. Even a single large batch will exceed
  this long before then.
- The hit count is `atomicAdd` on `unsigned int` and returned as `int`. Past
  2^31 hits the return is negative and the demo's retry loop
  (`if (found <= capacity) break;`) silently accepts a truncated result.

Both are **silent** wrong answers, which is the category that matters most here.
Fix with A6's named types, widening index/count to 64-bit.

### C2 [S]. `decompose_prime` fails silently, and the failure deletes solutions

`gaussian.cc`: `sqrt_minus_one` gives up after trying `n < 1000` and returns 0;
`decompose_prime` then returns `false`; `point_set` returns `{}`; the norm
contributes **zero points** and is skipped without a word.

The probability of a genuine miss is astronomically small (half of all residues
are non-residues, so 1000 failures in a row is ~2^-1000). That is not the point.
The point is that the failure mode is indistinguishable from "this norm has no
solutions", which is exactly the failure the paper's §6.8 is written to prevent:
*a bug in an exclusion search produces nothing, forever, and is invisible in the
output by construction.*

Same class: `point_set` returns `{}` when handed a non-split prime.

Both should throw. There is no performance argument for the silent path — this
runs once per prime, not per norm.

### C3 [S]. The wraparound false-positive path is completely untested

This is the one I most want fixed. Above `kExactKeyBound` (`L1 > 3.037e9`) the
kernel's `uint64` comparison is a *filter*, and the entire soundness argument
rests on `verify_exact` rejecting the false positives it lets through. That path
has no test, and no test can reach it by accident — a random false positive has
probability ~2^-58.

It is reachable by construction. Verified:

```
nu       = [2^33, 3, 2^32, 2^63 - 3]      (four distinct values)
exact k1 = 73786976294838206473
exact k2 = 85070591730234615828950163710522949641
equal exactly      : false
equal mod 2^64     : true   (both reduce to 9)
```

So the kernel emits this quadruple and `verify_exact` must reject it. That is a
five-line unit test on the CPU side and a five-line one on the GPU side, and it
pins down the only property that makes large-`B` results defensible.

Related: `kExactKeyBound` is defined but **never read anywhere**. Nothing warns
when a run crosses it.

### C4 [C]. `meet_batch` cannot validate its own inputs

It trusts `offsets[num_norms]` for the length of `nu`, which the caller never
passes. A malformed or non-monotonic `offsets` array is an out-of-bounds device
read with no diagnostic. Add an explicit `nu_count` parameter and validate that
offsets are non-decreasing and bounded by it. Cheap, and it is the kind of check
that pays for itself the first time a batching bug appears.

### C5 [C]. `decompose_prime` does not check primality

Handed a composite `n ≡ 1 (mod 4)` that happens to be a sum of two squares
(25 = 3²+4², say), it can return a perfectly well-formed decomposition of a
non-prime. The contract says "prime"; nothing enforces it, and the result is
silently wrong point sets. At minimum a debug assert.

### C6 [C]. `Norm::primes` aliases live DFS state

`for_each_admissible_norm` hands the callback a `const Norm&` whose `primes`
vector is the DFS's working stack — it is mutated as soon as the callback
returns. The demo happens to copy it; a caller that stores the reference gets
garbage. Undocumented. Either document loudly or hand out a `std::span` plus a
value type.

### C7 [P]. `for_each_admissible_norm` type-erases its callback

`std::function` in what will be a 2×10^10-iteration loop. Template the callback.

### C8 [P]. `split_primes_up_to` allocates a full sieve — hard wall at ~10^10

`std::vector<bool> composite(limit + 1)`. At `B = 10^12` the limit is
`B/65 ≈ 1.5e10`, so ~1.9 GB and minutes of sieving before the first norm is
emitted. `for_each_admissible_norm(1e12)` as written will not run. Segmentation
is already on the deferred list (largest-prime-outer); the *current* ceiling
should be documented in the header now so nobody discovers it the hard way.

### C9 [P]. Both meets are O(D⁴), with no bucketing

`meet_reference` and the kernel both walk 4-subsets: `C(D,4) × 3 ≈ D⁴/8`. At
`D = 36` (the `B = 10^8` maximum) that is 177k operations for one norm — fine.
At `D = 192` (the `B = 10^12` maximum) it is 1.6×10^8 for a *single norm*, which
is seconds. The algorithm doc's crossover to hashing at `D ≈ 8–10` is not
optional at scale. Known and deferred, but the number is worth having written
down.

### C10 [P]. Small transfer/allocation waste in `meet_batch`

`d_hits.CopyToHost` copies the whole `max_hits × 5` buffer even when one hit was
written; every call re-allocates all four device buffers. Subsumed by A2/A4 if
those are done properly.

### C11 [C]. Test gaps

- **No unit test that the CPU reference and the GPU agree on real norms.** The
  demo does this, but `ctest` does not, so a regression is only caught by
  running the demo by hand. A small-bound version belongs in the test suite.
- **The `bound/65` factor-base cutoff is load-bearing and untested.** I audited
  it externally by brute force: at `B = 10^6` the largest prime factor actually
  occurring in an admissible norm is **15377**, against a cutoff of
  `10^6/65 = 15384`. It is correct — and *exactly* tight, the extremal norm
  being `5 · 13 · 15377 = 999505`. A one-off-by-one in that reasoning silently
  deletes solutions, so it deserves an explicit test rather than a comment.
- Nothing tests `verify_exact` on values that need more than 64 bits (see C3).

### C12 [A]. The demo holds everything in memory

`enum_demo.cc` builds `work` for every norm before touching the GPU, so memory
is O(total points) and there is no batching. Fine at `B = 2×10^6`; it is the
next thing to break as the bound rises.

---

## Part C — consolidated TODO, in the order I would do it

**Tier 1 — silent-wrong-answer risks. These first.**

1. **C3** Test the wraparound/false-positive path, both CPU and GPU, using the
   constructed fixture above. Make `kExactKeyBound` actually load-bearing:
   assert or warn when a run crosses it.
2. **C2** Make `decompose_prime` and `point_set` throw instead of returning
   empty. No silent skips anywhere in the enumeration path.
3. **C1 + A6** Named types in a shared header; widen index/count to 64-bit;
   delete the `reinterpret_cast`s.
4. **C11** Add the CPU-vs-GPU unit test and the `bound/65` cutoff test.

**Tier 2 — architecture, best done as one pass.**

5. **A4/A5** Split `gpu.h` per app, with ownership hoisted into a session
   object rather than allocated per call; resolve the seam question (pimpl is
   my recommendation).
6. **A1** `check.cuh` → `check.cuh` + `device_buffer.cuh`.
7. **A7** `DeviceAppendBuffer<T>` with its own tests and an explicit overflow
   contract.
8. **A3** `#pragma once` sweep across all headers.
9. **C4** `nu_count` parameter and offset validation.
10. **C5, C6** Primality assert; document or remove the `Norm::primes` aliasing.

**Tier 3 — performance, gated on a measurement.**

11. **A2** Benchmark pageable vs. pinned vs. async *first*, then decide how much
    `HostBuffer` machinery to build. Producers must fill it directly.
12. **A8 + C9** `__restrict__` now; precomputed squares and `D`-bucketing
    together.
13. **C7, C8, C10, C12** Template the callback, segment the sieve, trim the
    copies, batch the demo.

**Open question worth settling before Tier 2.** Everything in Tier 2 and 3 is
shaped by whether the enumerator stays a naive hybrid or moves to
device-resident generation. The demo measured 0.25 ops/byte, and the algorithm
doc's §2.2 argues the fix is to generate norms inside the block and stop
transferring them. If we commit to that, a good deal of the transfer-path work
(A2 in particular) is optimizing a path we intend to delete. I would settle that
question first and let it set the scope of the refactor.
