# Code review — responses and TODO

**Status: Part A (style and convention) is done and verified; Parts B and C are
open.** Build is warning-free, all 8 test binaries pass, both demos produce
identical results to before the refactor.

Working file for the review of `code/native`. Part A answers the inline
`GILRO:` comments one by one; Part B is my own audit of the same code, which
turned up things the inline comments did not; Part C is the consolidated,
prioritized TODO.

Severity labels: **[S]** silent-wrong-answer risk (the failure mode that
actually threatens an exclusion claim), **[C]** correctness/robustness,
**[P]** performance, **[A]** architecture/hygiene.

Parts A, B and C below were written before any code changed, and are left as
written so the reasoning can be reviewed against the outcome; what actually
shipped is recorded at the end. Every claim with a number attached was checked.

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


---

## Round 1 outcome — what shipped

Everything in Part A except the two items explicitly deferred.

| Item | Disposition |
| --- | --- |
| A1 `DeviceBuffer` out of `check.cuh` | done — `core/device_buffer.cuh`; `check.cuh` keeps error handling |
| A2 `HostBuffer<T>` pinned | done — `core/host_buffer.cuh`, used by `meet.cu` for staging. Exposes `begin()/end()/operator[]` so a producer can fill pinned memory directly |
| A3 `#pragma once` | done — every header |
| A4/A5 split the monolith | done — `gpu.h` deleted; per-app directories, each with its own library and test binary |
| A6 named types | done — `core/types.h` (`nu_t`, `key_t`, `index_t`, `count_t`); every `reinterpret_cast` gone |
| A7 output buffer as a class | done — `core/device_append_buffer.cuh` with its own tests. `MeetHit` replaced the 5-int stride, so no manual index arithmetic survives anywhere |
| `__launch_bounds__` | done — on all three kernels |
| `PTE_CUDA_LAUNCH_KERNEL` | done — every launch goes through it |
| kernel naming | done — `E4MeetDirectKernel`, with the header comment spelling out both senses of "specific" |
| A8 amortize squarings | **deferred at your request** |

Two things rode along because they were inseparable from A6, and both are
strict improvements rather than behaviour changes:

- `index_t`/`count_t` widened to 64-bit, closing **C1**. Picking a width was
  unavoidable once the aliases existed, and leaving 32 bits in place would have
  entrenched a silent overflow in the signatures being rewritten.
- `meet_batch` gained a `nu_count` parameter and now validates that offsets are
  non-decreasing and in range, closing **C4**. Previously a malformed offsets
  array was an unchecked out-of-range device read.

### One regression worth watching

The demo's GPU stage went from 0.92 ms to 3.77 ms. That is `cudaHostAlloc`
being called per `meet_batch` invocation for the staging buffer — exactly the
"pinned allocation is expensive and must be pooled, never per call" caveat from
A2. It is invisible at the demo's size and would be ruinous at campaign scale.
The fix is the session object from A4/A5 owning the staging buffers across
calls; until then, treat the meet timing in the demo as not meaningful.

### Structural decision made along the way

Applying "each app in its own directory" to `src/gpu` exposed a latent hazard:
`cpu/bigint/bigint.h` and `gpu/bigint/bigint.h` share a basename, and with both
trees on the include path a bare `#include "bigint.h"` resolved by include-path
order. The include root is therefore now `src/`, and every header is reached by
its full path. `src/cpu` got the same per-app treatment for consistency.

### Still open

Parts B and C above, unchanged, minus C1 and C4. The three highest-value items
remain the silent-failure risks, none of which are style:

1. **C3** — the wraparound false-positive path is still untested, and
   `kExactKeyBound` is still never read anywhere.
2. **C2** — `decompose_prime` and `point_set` still fail silently.
3. **C11** — no CPU-vs-GPU unit test; the `bound/65` cutoff is still untested.

---
---

# Round 2

Responses to the second batch of inline `GILRO:` comments. Same convention as
Part A: agreement where I agree, pushback where I don't, and a number attached
to every claim that has one.

Two comments asked for a design rather than a fix ("let's discuss what should be
the proper separation here", "suggest a refactor together with `enumerate/`"), so
those get a proposal rather than a patch.

## Part D — responses

### D1. `.cuh` → `.h`, and `meet.cu` → `meet.cc` — **agreed, all six**

`check.cuh:11`, `device_buffer.cuh:11`, `host_buffer.cuh:32`, `meet.cu:16`,
`meet_impl.cuh:21`, and the separation comment at `device_append_buffer.cuh:25`.

You are right on every one of them, and the reason is worth writing down as a
rule because I got it wrong by default:

> **`.cu`/`.cuh` iff the file contains device code** — a `__device__`/`__global__`
> definition, `<<< >>>`, or a device intrinsic like `atomicAdd`. Everything
> else is `.cc`/`.h`, *including every call to the CUDA runtime API*.

`cudaMalloc`, `cudaMemcpy`, `cudaHostAlloc`, `cudaGetErrorString` are a C
library. They need `<cuda_runtime.h>` and `-lcudart`; they do not need nvcc. I
had been treating "mentions CUDA" as the test, which is why four host-only files
were `.cuh`.

The rule matters more than it looks, and this is the part I'd underweighted:
**one `__device__` function in a header promotes every translation unit that
includes it, transitively, to nvcc.** `DeviceAppendView::push` was doing exactly
that — it dragged `device_append_buffer` → `meet_impl` → `meet` into nvcc, which
is why `meet.cu` had to be a `.cu` even though it contains no device code at all.
That is the coupling your separation comment identifies, and it propagates.

What that bought once fixed, measured on this tree:

- `pte_gpu_meet` now has **one** nvcc translation unit (`meet_impl.cu`) instead
  of two.
- Three files moved to the host compiler — `meet.cc`, `meet_test.cc`,
  `host_buffer_test.cc` — and thereby **gained `-Wpedantic`**, which CUDA
  translation units in this project do not get (see `pte_warnings` in the root
  `CMakeLists.txt`: CXX gets `-Wall -Wextra -Wpedantic`, CUDA gets only
  `-Xcompiler=-Wall,-Wextra`).
- `meet_test.cc` being a `.cc` is now *mechanical* proof of the seam rule: a
  test that cannot contain device code cannot accidentally reach past `meet.h`.

Not done, and worth knowing: `bigint/bigint_kernel.cu` and
`matmul/matmul_kernel.cu` still mix host wrappers and kernels in one `.cu`. They
are scaffold kernels and the split is mechanical; I left them rather than widen
this round.

### D2. The view/kernel separation — **agreed, and there is a trap in the obvious version**

`device_append_buffer.cuh:25`. The split is: POD in the `.h`, the one device
operation in the `.cuh`.

```
core/device_append_buffer.h     DeviceAppendView<T> (data members only)
                                DeviceAppendBuffer<T> (the host owner)
core/device_append_buffer.cuh   __device__ push(view, value)
```

Host code needs to *name* and *construct* a view — `meet.cc` does
`d_hits.view()` and passes it to the launcher — and both are pure C++. Only
*writing through* it needs nvcc.

The trap, which is why `push` became a free function rather than staying a
member: the natural fix is to keep it a member and guard it with
`#ifdef __CUDACC__`. That gives the class **a different definition in host and
device translation units** — a textbook ODR violation, and one no compiler will
diagnose. A free function has neither problem: one type everywhere, only the
operation is device-side.

Your "could add `__host__`" on the same line: **not yet, and not this way.** A
`__host__ __device__` `push` would still call `atomicAdd`, which has no host
implementation, so the host half would need a `std::atomic` version operating on
host memory — a different object, not the same view. The real version of that
idea is a host type satisfying the same concept, so one kernel body can be
compiled for both. That is a genuine technique and worth doing when we have a
kernel body we want to share; it is not a one-line addition.

### D3. Scope of `types.h` — **agreed, and it should be narrower still**

`types.h:17`. Right — `nu_t` and `key_t` are meet concepts, not infrastructure.
Checking actual usage, the split is sharper than "move them to meet":

| type | used by | now lives in |
| --- | --- | --- |
| `index_t` | meet **and** the append buffer's capacity | `core/types.h` |
| `count_t` | the atomic counter only | `core/types.h` |
| `nu_t` | meet's public signature | `meet/meet.h` |
| `key_t` | **the kernel body, and nothing else** | `meet_impl.cu`, file-local |

The key type never crosses any API. It had no business being visible at all,
let alone project-wide.

`core/types.h` is also where I put the extension rule from D1, since that is
where a reader goes looking for conventions.

### D4. `key_t` → `trunc_key_t` — **agreed**

`types.h:25`. Better name, for the reason you give. It is now file-local in
`meet_impl.cu` with the danger spelled out at the declaration: above
`kExactKeyBound` a match on this type is a *candidate*, and the only thing
between it and a wrong answer is the host re-check.

(`key_mod64_t` would be more literally accurate — it is a reduction, not a
truncation of a computed value — but "trunc" carries the warning better, which
is the whole point of the rename.)

### D5. `MeetHit` sharing `pte::cpu::Quadruple` — **pushing back**

`meet.h:19`. (The type you mean is in `enumerate.h`, not `gaussian.h`.)

My objection is not the duplication — it is the **dependency edge**. Sharing the
type makes `pte_gpu_meet` link `pte_cpu_enumerate`, and:

- `pte_cpu` pulls in GMP and Boost. The GPU half stops being buildable or
  testable on its own.
- The CPU enumerator is the **oracle** the kernel is checked against. The value
  of `meet_test.cc` is that it links `pte_gpu_meet` *alone* and writes its own
  host loop, so the two implementations cannot share a bug. Pointing them at a
  common header weakens the fixture that the whole CPU-vs-GPU cross-check rests
  on.

They are also not the same type: `MeetHit` carries a norm index because the GPU
works in batches; `Quadruple` does not because the CPU works per norm.

If a third consumer ever wants one type, the right home is a neutral
`src/core/` that both halves depend on — not `cpu/`. For five fields I would not
pay for that yet. The conversion belongs in the app layer (`enum_demo.cc`),
which is allowed to depend on both.

### D6. Launch bounds, both arguments, as template parameters — **agreed, and the measurement changed my answer**

`meet_impl.cu:33`. Done: `E4MeetDirectKernel<kThreads, kMinBlocksPerSM>`.

I was going to leave the second argument at 1 and call that "unconstrained".
**That is wrong, and `-Xptxas -v` says so.** Measured on sm_89, 128 threads/block:

| `minBlocksPerSM` | registers | spilled | blocks/SM | occupancy |
| --- | --- | --- | --- | --- |
| *(absent)* | 56 | 0 | 9 | 75% |
| 1 | 66 | 0 | 7 | 58% |
| 6 | 62 | 0 | 8 | 67% |
| 9 | 56 | 0 | 9 | 75% |
| 12 | 40 | **68 B** | 12 | 100% |

`1` is not "no constraint" — it *promises only one resident block*, so nvcc
happily spends 66 registers and occupancy drops to 58%. `12` buys full occupancy
by spilling, which for a loop this ALU-heavy is the wrong trade.

Set to **9**: reproduces what nvcc picks unaided, and pinning it means a later
edit to the kernel body cannot silently drift to the 7-block configuration.
Caveat recorded at the constant: it is tuned for sm_89 specifically, and
occupancy here is a proxy, not a timing — the kernel is microseconds against a
transfer-bound demo, so these rows are not yet distinguishable in wall clock.

### D7. Fixed-width types everywhere, GMP for validation — **agreed on the principle; the boundary is not where you put it**

`enumerate.cc:20`, `gaussian.cc:11`, `gaussian.cc:18`.

The principle is right and it is the single most important thing in this review:
**an exclusion claim dies from silent overflow**, so the project needs a visible
split between a fast path that is allowed to be a *filter* and an exact path
that is the *arbiter*, with a stated soundness relation between them. The meet
already has that (uint64 wraparound + `__int128` re-check); your comment is
asking to make it a policy instead of an accident. Agreed, and it should live in
one header that states the bound for every quantity.

Where I disagree is the claim that GMP is what the exact path needs. The
arithmetic, at `B = 10^12`:

| quantity | bound | bits | exact in |
| --- | --- | --- | --- |
| Gauss coordinates `m, n` | `sqrt(B)` = 10^6 | 20 | **int32** |
| `L1`, `nu = m²−n²` | `B` | 40 | int64 |
| pair key `nu_a² + nu_b²` | `2B²` | **81** | **`__int128`** |

The exactness thresholds for the pair key:

```
uint64  (2^64):  L1 < 3_037_000_499        <- kExactKeyBound, today's filter boundary
__int128 (2^127): L1 < 9_223_372_036_854_775_808
```

So `__int128` is exact for the pair key up to `L1 < 9.2×10^18` — **nine million
times past `B = 10^12`**. There is no plausible campaign in which the meet or
the point-set layer needs GMP, and putting it there would cost a heap allocation
per operation in the one loop that must not have one.

GMP *is* required, and we will need it — one layer further out. Verifying a
recovered ideal symmetric PTE of size 16 means power sums to `k = 15` of 16
integers of size `~sqrt(2·L1) ≈ 1.4×10^6`: **~311 bits**. Well past `__int128`.
A size-32 identity would need ~638.

So the boundary I'd propose is:

```
enumeration + meet          int32/int64 fixed width, wraparound where documented
candidate confirmation      __int128, exact to L1 < 9.2e18   (verify_exact today)
certificate / PTE identity  GMP                              (not yet written)
```

which is a *three*-way split, not two. The thing your comment is really pointing
at — "it should be clear" — I'd fix with a single `precision.h` naming each
quantity, its bound, and `static_assert`s tying the types to `B`. That kills the
`using u64` buried in a `.cc` (`gaussian.cc:11`) at the same time, which is fair
criticism: that alias is load-bearing and invisible.

### D8. Templating `Gauss` on the integer type — **agreed, but the payoff runs the other way**

`gaussian.h:20`. Worth doing. The motivation you give is future width; the
near-term win is the opposite direction. Coordinates are `< 10^6` — 20 bits —
so `Gauss<int32_t>` **halves the bytes** of every point. In a pipeline the demo
measures at **0.25 ops/byte**, transfer-bound by a factor of ~80, that is not a
micro-optimization.

One trap to get right, and it is exactly the class of bug the exact/fast split
exists to prevent: `norm()` and `nu()` of a `Gauss<int32_t>` **do not fit in
int32** (they reach `B` = 10^12). The template needs a separate, wider result
type or it silently overflows the moment anyone instantiates it narrow:

```cpp
template <typename T> struct Gauss {
  using wide = std::conditional_t<sizeof(T) <= 4, std::int64_t, __int128>;
  T re{}, im{};
  wide norm() const { return wide{re} * re + wide{im} * im; }
  wide nu()   const { return wide{re} * re - wide{im} * im; }
};
```

### D9. Member functions and operator overloading — **agreed**

`gaussian.h:33/36/39/69`. `canonize()`, `conj()`, `norm()`, `nu()` as members and
`operator*` for `mul` — all better, and they fall out of D8 anyway.

One exception: **`mul_and_conj_mul` stays a free function.** It returns two
results from four multiplies (both `z·q` and `z·q̄`), which is the fused trick at
the centre of point-set construction. There is no operator for that shape, and
splitting it into two operators throws away the saving it exists for.

### D10. `point_set` recomputes decompositions; wants an `Enumerator` — **agreed, and it is a real cost, not just shape**

`gaussian.cc:101`, with `enumerate.cc:35` and `:70`. This is the best catch in
the round. Measured:

| bound | admissible norms | `decompose_prime` calls | distinct primes | redundancy |
| --- | --- | --- | --- | --- |
| 10^6 | 11,083 | 36,675 | 891 | **41×** |
| 10^7 | 131,437 | 436,126 | 7,072 | **62×** |

And directly, at the demo's `B = 2×10^6` (78,047 calls):

```
point_set total          8.12 ms
  decompose_prime only   4.90 ms   = 60.4% of it
```

**Three fifths of point-set construction is recomputing the same few thousand
Cornacchia descents.** Point-set construction is in turn the largest CPU cost in
the demo (11.8 ms of ~19 ms).

Two levels of fix, and the second is the one you're describing:

1. **Memoize.** The factor base is already materialized; carry a parallel
   `std::vector<Gauss>` of decompositions. Cuts that 60% to nothing. Half a day.
2. **The `Enumerator` class.** The deeper waste is that each norm rebuilds its
   *entire* point set from scratch, when a DFS node's point set is its parent's
   doubled by one prime. An object owning the factor base, the decomposition
   table, and a stack of partial point sets along the DFS spine makes the
   incremental structure available. This is exactly §A2 (largest-prime-outer) in
   `drafts/gpu_algorithms.md` — so your instinct here and the written plan agree,
   which is reassuring for both.

I'd sequence memoization first (cheap, immediate, no design risk) and fold the
`Enumerator` into the CPU-implementation round, since it is the natural owner of
the batching that C12 wants anyway.

### D11. `std::function` → template — **agreed; honest about the size**

`enumerate.cc:70`. This is C7 from Part B, so we agree. Templating the callback
also lets the DFS body inline, which `std::function` blocks entirely.

Sizing it honestly, though: the indirect call is ~2 ns against ~700 ns of
`point_set` work per norm, so **this is not measurable today** — `point_set`
(D10) is where the time is. It becomes real at the 2×10^10 iterations the
campaign implies. Worth doing, worth doing after D10.

### D12. `verify_exact` "should live under test" — **disagreeing, and this one matters**

`enumerate.cc:149`. `verify_exact` is not a test helper. Above `kExactKeyBound`
(`L1 > 3.037×10^9`) the kernel compares keys modulo 2^64 and therefore emits
**candidates**; `verify_exact` is the step that turns candidates into results.
Move it to test code and the production path silently accepts false positives —
which is precisely **C3**, the highest-severity open item.

Your instinct is still telling you something true: it is *named* like a test
helper, and today it is only ever called by the demo, so nothing forces it. Both
of those are the actual defect. The fix is the other direction:

- rename it to something that carries its role — `confirm_candidate`, or
  `is_true_collision`;
- make it **unskippable**, by having the enumerator return confirmed hits only,
  so an unverified candidate cannot escape the API.

### D13. `decompose_prime` should return `Gauss` and throw — **agreed**, and it is the same bug as the magic number

`gaussian.cc:56` and `:66`, plus **C5** from Part B. These three are one defect.

`sqrt_minus_one` trials `n < 1000` for a non-residue. Half of all residues are
non-residues, so for a genuine prime `p ≡ 1 (mod 4)` the chance of 1000
consecutive failures is ~2^-1000. It cannot happen. Which means:

> **if that loop ever exhausts, `p` is not prime.**

The magic number and the missing primality check are the same thing — the loop
*is* an (accidental, unreadable) primality test whose failure is currently
swallowed into `return 0` → `return false` → `point_set` returns `{}` → the norm
contributes zero points and is skipped in silence. In an exclusion search that
is the worst available failure mode.

So: `Gauss decompose_prime(u64 p)`, throwing `std::invalid_argument` on a
non-split-prime input, with the bound named (`kMaxNonResidueTrials`) and the
2^-1000 argument at the declaration. Not `std::optional` — the factor base
contains only split primes by construction, so a failure here is a programming
error at the call site, not an expected outcome.

### D14. `isqrt` and GMP — **pushing back on the reason, agreeing with the instinct**

`gaussian.cc:26`. Two things in the comment, and they point opposite ways.

"Reinventing the wheel **in a very inefficient manner**" — this one I think is
backwards. The current `isqrt` seeds from `std::sqrt(double)` and corrects,
which is one or two iterations for any 64-bit input: tens of nanoseconds.
`mpz_sqrt` on a value this small costs an `mpz` init and a library call, which
is several times *more*. Swapping it in would slow down a function that D10 just
showed is 60% of point-set time.

But "don't trust a hand-rolled `isqrt`" is right, and has a track record here —
an earlier version of this exact function in this project had a malformed
correction loop. The resolution that gets both: **keep the fixed-width version,
and test it against `mpz_sqrt`.** GMP is already a dependency (`pkg_check_modules
... gmp gmpxx`), so a property test over random `u64` plus the boundaries
(perfect squares, `n²−1`, `2^64−1`) costs nothing and gives us the wheel's
authority without its call overhead. That belongs in the same test round as C11.

### D15. `std::span` instead of raw pointers — **agreed, and it deleted four tests**

`meet.h:66`. Right, and it pays for itself immediately. Before and after:

```cpp
index_t meet_batch(const nu_t *nu, index_t nu_count, const index_t *offsets,
                   index_t num_norms, MeetHit *hits, index_t max_hits);

index_t meet_batch(std::span<const nu_t> nu, std::span<const index_t> offsets,
                   std::span<MeetHit> hits);
```

Six parameters to three. Each of the three that vanished was a *count the caller
could get wrong independently of the data it described* — and note that **C4,
the out-of-bounds device read from round 1, was exactly that bug**: `nu_count`
had to be added last round precisely because `meet_batch` could not otherwise
know how long `nu` was. With a span the question cannot be asked. `num_norms`
likewise becomes `offsets.size() - 1`, enforced by the type rather than by a
comment.

The visible measure: `RejectsBadArguments` lost four of its six cases —
negative `num_norms`, negative `max_hits`, null `nu`, null `hits` — because they
can no longer be *written*. What survives is the one class a span cannot cover:
`offsets` is a *description of* `nu`, so only a cross-check between the two
catches a malformed one. That test is renamed `RejectsMalformedOffsets` and
keeps the three real cases (past the end, non-monotonic, negative).

Two notes on the edges:

- **This required `CMAKE_CUDA_STANDARD` 17 → 20.** `meet.h` is included by
  `meet_impl.cu`, and `std::span` is C++20. Verified nvcc 12.9 compiles C++20
  with `<span>` before making the change; host code was already C++20.
- **Raw pointers stay at the device contract** (`meet_impl.h`,
  `LaunchE4MeetDirect`). Those address device memory, and `std::span`'s members
  are not `__device__`, so a span cannot be indexed inside a kernel. The honest
  boundary is: std types at the host API, raw pointers at the device seam.

"Lvalue refs" I read as the same request; for array arguments `std::span` is the
form that carries it, and there is no scalar out-parameter left to take by
reference.

---

## Round 2 outcome — what shipped

Everything in D1–D6 and D15, which is the whole GPU-side batch. D7–D14 are CPU-side and
are a redesign (`Enumerator`, templated `Gauss`, `precision.h`); those are
written up above as a proposal and **not implemented** — they want your call
first, and they land naturally in the CPU-implementation round.

| Item | Disposition |
| --- | --- |
| D1 extensions | done — `check.h`, `device_buffer.h`, `host_buffer.h`, `meet_impl.h`, `meet.cc`, `meet_test.cc`, `host_buffer_test.cc`; rule documented in `core/types.h` |
| D1 launch macro | `PTE_CUDA_LAUNCH_KERNEL` split into `core/launch.cuh` — `<<< >>>` is the one thing in `check.cuh` that actually needed nvcc |
| D2 view/kernel split | **superseded — see D16.** `push` is a `__device__` member under `#ifdef __CUDACC__`; the `.cuh` is gone |
| D2 `__host__` on push | not done — would need a host `atomicAdd`; see above |
| D16 `LaunchE4MeetDirect` takes `DeviceBuffer&` | done — raw device pointers now appear exactly once, in the kernel argument list |
| D3 type scope | done — `nu_t` → `meet.h`, `key_t` → file-local in `meet_impl.cu`, `core/types.h` keeps only `index_t`/`count_t` |
| D4 `trunc_key_t` | done |
| D5 share `Quadruple` | **not done — pushing back**, see above |
| D6 launch bounds | done — both as template parameters, `minBlocksPerSM = 9` from measurement |
| D15 `std::span` API | done — `meet_batch` 6 params → 3; `CMAKE_CUDA_STANDARD` 17 → 20 |
| D7–D14 | proposals above; awaiting your call |

Verification: warning-free build, **8/8 test binaries pass**, both demos
unchanged (same two solutions at `L1 = 67405` and `600445`, same one imprimitive
rejection, CPU and GPU sets identical).

Also closed, carried over from round 1: **`compute-sanitizer --tool memcheck` is
clean on all four GPU test binaries.**

Still not done from round 1: the pageable-vs-pinned bandwidth demonstration in
`hybrid_demo.cc`, and the per-call `cudaHostAlloc` regression (0.92 ms → ~4 ms)
that needs the session object to fix.

## Part E — updated TODO

**Tier 1 — silent-wrong-answer risks. Unchanged, and still first.**

1. **C3** Test the wraparound false-positive path with the constructed fixture
   (`nu = [2^33, 3, 2^32, 2^63−3]`); make `kExactKeyBound` load-bearing.
2. **C2 + D12 + D13** One change, not three: `decompose_prime` returns `Gauss`
   and throws; the trial bound gets a name and doubles as the primality signal
   (**C5**); `verify_exact` is renamed and made unskippable.
3. **C11 + D14** CPU-vs-GPU unit test; the `bound/65` cutoff test; `isqrt`
   against `mpz_sqrt`.

**Tier 2 — the CPU implementation round.**

4. **D7** `precision.h`: name every quantity, its bound, and the three-way
   exact/filter/certificate split, with `static_assert`s.
5. **D10** Memoize `decompose_prime` (cheap, 60% of point-set time), then the
   `Enumerator` owning factor base + decomposition table + DFS point-set stack.
6. **D8 + D9** Template `Gauss`, with the wider result type; members and
   `operator*`; `mul_and_conj_mul` stays free.
7. **D11 / C7** Template the callback once `Enumerator` exists.
8. **C6** `Norm::primes` aliasing — document or hand out a value type.

**Tier 3 — performance, gated on the architecture decision.**

9. Session object owning staging buffers (fixes the per-call `cudaHostAlloc`).
10. **A2** Pageable-vs-pinned benchmark, then decide how much staging machinery
    to keep.
11. **A8 + C9** Precomputed squares and D-bucketing together.
12. **C8, C10, C12** Segment the sieve, trim the copies, batch the demo.

**Still the open question** from round 1, now more pointed: naive hybrid or
device-resident generation. D8 (narrow `Gauss` halving transfer bytes) is worth
real money under the first and nearly nothing under the second.


---

## Round 3 — two corrections and a closing note

### D16. `push` as an `#ifdef __CUDACC__` member — **your call, and it is better than mine**

You overruled the free function. Implemented as asked, and on reflection the
argument I made in D2 was weighing the wrong thing.

What I got right: the guard does give the class different definitions in host
and device translation units, and that is formally an ODR violation.

What I underweighted: **it is bounded and benign here, and it deletes a file.**
The members are identical, so layout and the by-value kernel argument are
unaffected; `push` is `inline __device__`, so it is only ever emitted into device
code; there is no virtual function or static member for the two definitions to
disagree about. It is also the idiom Thrust and CUB are built on. Against that,
the free-function version cost an extra header *and* an unnatural call site —
and your original complaint was precisely that the structure was forcing `.cuh`
files where `.h` would do. The guard serves that goal better than my split did:

```
before:  device_append_buffer.h  +  device_append_buffer.cuh   (push(out, v))
after:   device_append_buffer.h                                (out.push(v))
```

One invariant now written at the declaration, because it is what makes this
safe: **anything added inside the guard must not touch layout.**

### D17. `LaunchE4MeetDirect` should take `DeviceBuffer`, not raw pointers — **agreed**

`meet_impl.cu:78`. Right, and it is the same argument as D15 one layer down:
`DeviceBuffer<T>` carries its own count, so `num_norms` stops being a parameter
a caller can contradict and becomes `d_offsets.count() - 1`.

```cpp
void LaunchE4MeetDirect(const DeviceBuffer<nu_t> &d_nu,
                        const DeviceBuffer<index_t> &d_offsets,
                        DeviceAppendView<MeetHit> out);
```

Four parameters to three, and the invariant "offsets has one more entry than
there are norms" is now expressed once, in the type, instead of being restated
by every caller. Raw device pointers survive in exactly one place in the app —
the kernel's argument list, which cannot take anything else.

This retires the "raw pointers at the device seam" caveat I attached to D15;
that was true of `std::span` specifically, not of the seam.

### GPU part: closed

`src/gpu` is finalized. Final state: **three `.cu` files** in the whole tree
(`device.cu`, `bigint_kernel.cu`/`matmul_kernel.cu`/`meet_impl.cu`), one
`.cuh` (`launch.cuh`), everything else ordinary C++.

Verified: warning-free clean build, 8/8 tests, `cpu-only` preset 4/4, both demos
unchanged (same two solutions at `L1 = 67405` and `600445`), 56 registers with 0
spills, compute-sanitizer memcheck clean.

Not carried further, and still true: `bigint_kernel.cu` and `matmul_kernel.cu`
still mix host wrappers with kernels in one `.cu`. They are scaffold code and
the split is mechanical.


---
---

# CPU implementation round

D7–D14 built. What follows records what shipped, the two places measurement
contradicted a claim I had made, and what is left.

## What shipped

| Item | Where |
| --- | --- |
| **D7** precision policy | `cpu/core/precision.h` — the three-way split, named types, `static_assert`s tying them to `kMaxNorm`; `kExactKeyBound` moved here from `enumerate.h` |
| **D8/D9** `GaussT<T>` | templated on coordinate type, with `wide` result type; `conj/canonize/norm/nu/divisible_by` as members, `operator*`; `mul_and_conj_mul` stays a free-standing two-output method |
| **D10** `Enumerator` | owns the factor base and its decompositions; `point_set(Norm)` uses them |
| **D11** templated callback | `for_each_norm(F&&)`, `std::function` gone |
| **D12** confirmation | `verify_exact` → `confirm_candidate`, documented as a *production* step, not a test helper |
| **D13/C2/C5** `decompose_prime` | returns `Gauss`, throws `NotASplitPrime`; deterministic Miller-Rabin |
| **D14/C11** tests | `isqrt` and `is_prime` against GMP; the `bound/65` cutoff by brute force; the wraparound fixture; `Norm` copy contract |
| **C11** cross-check | `enum_demo` registered as ctest `enum_demo_cross_check` |

Test binaries: 8 → 10. Both presets warning-free; `dev` 10/10, `cpu-only` 5/5;
demo unchanged (`L1 = 67405`, `600445`).

## Two things measurement corrected

### The trial-division loop is not a primality check — my C5 argument was wrong

I claimed in D13 that exhausting `kMaxNonResidueTrials` means `p` is composite,
and that this made the magic number and the missing primality check "one bug".
The first half is true; the conclusion is not. I wrote the test before the code
and it failed on **85, 125, 221**.

Those are composites *and* products of split primes, so they are sums of two
squares: `85 = 9² + 2²`, `221 = 14² + 5²`. A non-residue exists, the descent
closes, and the old code returned a well-formed decomposition of a non-prime.
And these are exactly the composites a caller is likely to pass by accident — a
partially-factored norm, say.

So C5 needed a real primality test, not a clever reading of an existing loop.
`is_prime` is now deterministic Miller-Rabin: `{2,3,5,7}` below 3,215,031,751
and `{2,3,5,7,11,13}` below 3,474,749,660,383, which covers the whole factor
base at `B = 10^12` (`1.5e10`). Exact, not probabilistic.

Cost, measured at `B = 2×10^6`: the whole `Enumerator` constructor — sieve,
1,649 Miller-Rabin tests, 1,649 Cornacchia descents — is **0.85 ms**. It is free
at this scale. It will not be at `B = 10^12`, where the factor base is ~3×10^8
primes; but the sieve there is a 1.9 GB allocation (**C8**) and hits the wall
first, so segmentation is the blocking item either way.

### D10 delivered, and the number is better than the 60% I predicted

| | before | after |
| --- | --- | --- |
| `decompose_prime` calls | 78,047 | **1,649** |
| point-set construction | 8.12 ms | 3.63 ms |
| + factor base and decompositions | (in the above) | 0.85 ms |
| **total** | **8.12 ms** | **4.48 ms** |

**1.8× on the CPU's largest cost**, at `B = 2×10^6` over 23,576 norms. In the
demo, point-set time went 11.8 ms → 7.6 ms across five runs.

## What I did not do

- **`GaussT<std::int32_t>` is not instantiated anywhere.** The template is
  shaped for it and `wide` handles the result-type trap, but `coord_t` stays
  `int64_t`. Narrowing pays when point *coordinates* cross the bus; today the
  meet moves `nu` values, which need 40 bits. Flipping it is a one-line change
  in `precision.h` when device-resident generation lands.
- **`Norm::primes` still aliases live DFS state** (C6). A value type means a
  heap allocation per norm in a loop that will run ~2×10^10 times. The contract
  is now documented at `for_each_norm` and there is a test that a *copied* Norm
  survives, which is what a caller actually needs.
- **A2 largest-prime-outer** — carrying point sets down the DFS spine rather
  than rebuilding per norm. This is the next real speedup and `Enumerator` is
  now the object that will own it.
- **Layer 3 (GMP certificates)** — nothing yet verifies a recovered PTE
  identity. That is where GMP finally earns its place; see the 311-bit row in
  `precision.h`.

## Open, in order

1. **C8** Segment the sieve — the hard wall at `B ≳ 10^10`.
2. **A2** Incremental point sets along the DFS spine.
3. Session object owning GPU staging buffers (the per-call `cudaHostAlloc`
   regression, still unfixed).
4. **A8 + C9** Precomputed squares and D-bucketing.
5. **C12** Batch the demo; it still holds every norm in memory.


---

## Bug found at scale: the least-non-residue cap

Surfaced while measuring whether the factor base is a "wall" (it is not — see
`drafts/implementation_design.md` §3.1). Building the $B = 10^{10}$ factor base
aborted on:

```
decompose_prime(48473881): no quadratic non-residue found in 64 trials,
                           so p is composite
```

**48,473,881 is prime.** Its least quadratic non-residue is 67.

This was my regression, introduced in D13. The original code searched
`n = 2..999`; "naming the magic number" I replaced it with
`kMaxNonResidueTrials = 64`, justified by: *half of all residues are
non-residues, so 64 consecutive failures has probability 2^-64.*

**That argument is wrong.** The trials are not random samples — they are the
smallest integers 2, 3, 4, …, and the least quadratic non-residue `n(p)` is a
studied quantity that is routinely far above 2. Measured by exhaustive scan of
every split prime in each factor base:

| $P_{\max}$ | split primes | max $n(p)$ | attained at |
| --- | --- | --- | --- |
| $1.54\times10^6$ | 58,433 | 41 | 1,083,289 |
| $1.54\times10^7$ | 497,241 | 53 | 9,257,329 |
| $1.54\times10^8$ | 4,327,911 | **67** | **48,473,881** |

Roughly +13 per decade, consistent with $n(p) \sim \log^2 p$. So 48,473,881 is
the record-holder for the entire $B = 10^{10}$ factor base, and the cap of 64
landed just under it.

**Fix.** The bound is now a runaway guard (`1 << 16`), not a probabilistic
threshold, with the reasoning corrected at the declaration:

- termination is *guaranteed* because `is_prime` has already run — for a prime
  there are exactly $(p-1)/2$ non-residues, so the loop must find one;
- under GRH, $n(p) < 2\ln^2 p$ (Ankeny) ≈ 1,100 at $P_{\max} = 1.5\times10^{10}$,
  so the guard carries ~60× headroom over the strongest conjectural bound;
- exhausting it means a bug here or a counterexample to GRH, so it throws.

Regression test added covering all three record-holders.

### What this says about the design

The failure was **loud**, and that is the whole point of D13. Before that
change, `decompose_prime` returned `false`, `point_set` returned `{}`, and the
norm contributed zero points — so this bug would have silently deleted every
admissible norm divisible by 48,473,881, and every campaign above $B \approx
3\times10^9$ would have been quietly incomplete. Nothing in the output would
have shown it.

That is exactly the failure mode the exclusion claim cannot survive, and the
throw converted it into a crash on the first run at scale.
