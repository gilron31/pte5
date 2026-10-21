# `code/native` — the compiled side of the E₅ work

Four things live here:

1. **The E₄ enumeration** — [`src/cpu/gaussian/`](src/cpu/gaussian/),
   [`src/cpu/enumerate/`](src/cpu/enumerate/), [`src/gpu/meet/`](src/gpu/meet/)
   and [`apps/enum_demo.cc`](apps/enum_demo.cc). This is search **S-A** of
   [`drafts/gpu_algorithms.md`](../../drafts/gpu_algorithms.md), in its first
   working form — see the section below.
2. **The champion meet** — [`src/core/`](src/core/),
   [`src/cpu/champion/`](src/cpu/champion/),
   [`src/gpu/champion/`](src/gpu/champion/) and
   [`apps/champion_demo.cc`](apps/champion_demo.cc). This is search **S-E** of
   [`drafts/gpu_implementation_plan.md`](../../drafts/gpu_implementation_plan.md),
   milestones G3 and G4, and it is the one search where the GPU earns its keep.
3. **The local model** — [`src/cpu/tower/`](src/cpu/tower/) and
   [`apps/tower_demo.cc`](apps/tower_demo.cc). The depth-*d* tower recursion over
   ℤ/pᵏ with a cap on node valuations: milestone **G1** of
   [`drafts/local_model_plan.md`](../../drafts/local_model_plan.md), which is
   where σ_p, the sieve wheels and any nonexistence certificate come from. Pure
   CPU, and it links nothing at all — see the section below.
4. **The original scaffold** — big-integer add/multiply (GMP vs. a naive
   schoolbook kernel) and matrix multiply (plain loops vs. a 16×16 tiled
   kernel), kept because it proves the plumbing end to end:

| Capability | Where |
|---|---|
| C++ library against Boost and GMP | [`src/cpu/`](src/cpu/) → `pte_cpu` |
| CUDA kernels running on the local GPU | [`src/gpu/`](src/gpu/) → `pte_gpu` |
| CPU processing that offloads to the GPU | [`apps/hybrid_demo.cc`](apps/hybrid_demo.cc) |
| CPU and GPU unit tested **independently** | `*_test.cc` / `*_test.cu`, CTest labels `cpu` / `gpu` |
| Python bindings to both paths | [`src/bindings/`](src/bindings/) → `pte_native` |

## Quick start

```sh
cd code/native

cmake --preset dev            # configure (fetches GoogleTest on first run)
cmake --build --preset dev    # build everything

ctest --preset dev -L cpu     # CPU unit tests
ctest --preset dev -L gpu     # GPU unit tests
ctest --preset dev            # both

./build/dev/apps/enum_demo       # the E₄ enumeration demo
./build/dev/apps/champion_demo   # S-E, the champion meet (arg: prime count)
./build/dev/apps/tower_demo      # G1, the local model (args: max p^k, depth, p)
./build/dev/apps/hybrid_demo     # the original scaffold demo

pip install -e .              # build + install the Python module
pytest                        # Python binding tests
```

## Read this first: your `PATH` has the wrong `nvcc`

This machine has **two** CUDA toolkits:

- `/usr/bin/nvcc` — CUDA **12.0**, from Ubuntu's `nvidia-cuda-toolkit` package.
  Its `crt/host_config.h` hard-errors on GCC > 12, so it **rejects the default
  `g++-13`**.
- `/usr/local/cuda` → `/usr/local/cuda-12.9` — CUDA **12.9**, from NVIDIA's
  `cuda-toolkit-12-9`. Accepts GCC ≤ 14. This is the one you want.

Nothing on the system puts `/usr/local/cuda/bin` on `PATH`, so a bare `nvcc`
resolves to the 12.0 one and any naive build fails with a confusing
"unsupported GNU version" error.

**This project does not depend on your shell being right** —
[`CMakePresets.json`](CMakePresets.json) and [`pyproject.toml`](pyproject.toml)
both pin `/usr/local/cuda/bin/nvcc` and `/usr/bin/g++-13` explicitly. But you
will want the shell fixed for everything else. Add to `~/.bashrc`:

```sh
# Go through the /usr/local/cuda symlink, which update-alternatives manages, so
# switching toolkits later is just: sudo update-alternatives --config cuda
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
```

Removing the 12.0 package is *not* necessary once `PATH` is right, and is mildly
risky: `apt` then reports ~60 packages as orphaned, including `g++-12`,
`libstdc++-12-dev`, and the intel-oneapi MKL runtime. If you do remove it, do
**not** run `apt autoremove` afterwards.

## The E₄ enumeration

`./build/dev/apps/enum_demo` runs search S-A end to end over admissible
$L_1 \le 2\times10^6$ and exits non-zero if anything disagrees.

**What it computes.** Work at the *norm-$L_1$ level*: for each admissible $L_1$
(a product of split primes with at least four representations as a sum of two
squares) build the point set $U(L_1) = \{z : N(z) = L_1\}$ modulo units and
conjugation, take $\nu(z) = \mathrm{Re}(z^2) = m^2 - n^2$ for each point, and
look for four points splitting as $\{z_1,z_2\}\mid\{z_3,z_4\}$ with

```
nu1^2 + nu2^2 == nu3^2 + nu4^2
```

That single condition *is* the E₄ problem. A root pair is recovered from a point
by the twist $z \mapsto (1+i)z$, giving $(r,s) = (m-n, m+n)$ with
$r^2+s^2 = 2L_1$.

Three facts from the algorithm doc shape the code:

- **$\nu$ is injective** on a point set, since $m^2+n^2 = L_1$ and $m^2-n^2 = \nu$
  determine $(m,n)$. So two distinct pairs with equal key can never share a
  point — there is **no disjointness test** anywhere in the meet.
- **The key needs more than 64 bits** in general ($\nu_1^2+\nu_2^2 \le 2L_1^2$).
  The kernel compares modulo $2^{64}$ by unsigned wraparound, which is exact
  below `pte::cpu::kExactKeyBound` and a *sound filter* above it: equal keys are
  always congruent, so there are **no false negatives**, only false positives
  that `verify_exact` rejects on the host. That asymmetry is what makes an
  exclusion claim defensible.
- **Repeated primes need care.** For $L_1 = p^bM$ a point is divisible by $p$
  iff its $p$-exponent $c$ satisfies $1 \le c \le b-1$, so an imprimitive
  fraction $(b{-}1)/(b{+}1)$ of the point set is exactly $p\cdot U(M)$. A
  quadruple with **all four** points imprimitive is $p^2\times$ a solution at
  $M$ — a duplicate, rejected by `is_primitive`. A quadruple with only *some*
  imprimitive points is genuinely primitive and is kept: that is the case where
  one pair uses the rational factor $p = \mathfrak{p}\bar{\mathfrak{p}}$ while
  another uses $\mathfrak{p}^2$, and discarding it would silently break
  exhaustiveness.

**Why the demo is slower than it looks like it should be.** It is deliberately
the *naive hybrid* split: the CPU generates norms and point sets, the GPU only
does the meet. That moves ~8 bytes per point to do a handful of integer
operations — the demo prints the resulting arithmetic intensity (~0.25 ops/byte),
which is transfer-bound by a wide margin. Generating the norms on the device is
what fixes it; see Deferred below.

**Correctness.** The shipped fixture is the CPU-vs-GPU cross-check in the demo:
`meet_reference` and the kernel must produce identical solution sets. The unit
tests add a cross-*language* check — the admissible-norm counts must match the
independent numpy sieve in
[`drafts/verify_gpu_algos.py`](../../drafts/verify_gpu_algos.py) exactly
(11,083 norms / 99,354 pairs at $10^6$; 131,437 / 1,339,509 at $10^7$) — plus
the paper's exemplar norm $13\cdot89\cdot173\cdot233\cdot709$, which must yield
its published solution.

### Deferred

Recorded so the next revisions have a defined starting point; none of this is
built yet.

*Fixtures.* Bremner family self-test (drive `bremner_F1_point` /
`bremner_F2_point` from [`code/old/utils.py`](../old/utils.py) and re-find every
member **including non-squarefree $L_1(k)$** — the test that would expose a
coverage hole); synthetic solution injection; a standalone non-squarefree
primitivity suite.

*Surface.* CLI flags for the bound and for CPU-only/GPU-only/both, so the demo
doubles as the pilot tool; Python bindings through `pte_native`, enabling
cross-checks against `code/old/enumeration_v1.py`; catalog output (CSV/JSON)
carrying the recovered towers so results accumulate.

*Algorithmic revisions*, in the order they matter:

1. **Largest-prime-outer generation** — emit each norm once as $n = c\,p^e$,
   amortizing Cornacchia to once per prime instead of once per norm.
2. **$D$-bucketing** — register-only kernels for $D = 4, 6, 8$ (~96% of norms);
   the bucket is known before generation since $D(cq) = 2D(c)$.
3. **Device-resident generation** — one block per cofactor subtree with the
   ancestor point sets in shared memory. This is what fixes the intensity
   problem the demo demonstrates.
4. ~~Champion-regime kernel~~ — **built**; see the next section. The
   $\mathbb{F}_p$ tower sieve, the top-collision and bottom-lift joins, and the
   A4′ product-form investigation are still open, and the tower sieve is blocked
   on a mathematics question rather than on code (plan §2.2).

## The champion meet (S-E)

`./build/dev/apps/champion_demo [num_primes]` runs search S-E end to end and
exits non-zero if anything disagrees.

**What it computes.** One hand-chosen norm with an enormous point set — $k$
distinct split primes give $D = 2^{k-1}$, so twenty of them give $D = 2^{19}$ at
$L_1 \approx 2.5\times10^{36}$ — and every collision *inside* it. Same condition
as S-A, `nu_a^2 + nu_b^2 == nu_c^2 + nu_d^2`; only the scale differs, and the
scale changes everything:

- $\binom{D}{2} = 1.374\times10^{11}$ pair keys, which at 16 bytes each is 2.2 TB
  and cannot be materialized. They are cut into **value windows**: with `sq[]`
  sorted, the pairs whose key falls in `[lo, hi)` are a contiguous run of `j`
  for each `i`, so each pair is enumerated exactly once instead of once per
  window pass. That is a 1000× correction over the hash-windowing both earlier
  drafts assumed.
- $\nu$ is 122 bits and the key ~244, so neither `int64` nor `__int128` is wide
  enough. [`src/core/uint.h`](src/core/uint.h) is the fixed-width type, shared
  by the host and the device.
- Point sets at that width are multi-limb Gaussian arithmetic
  ([`champion_points.cc`](src/cpu/champion/champion_points.cc)), not the `int32`
  path — whose coordinates stop at 20 bits.

**Two implementations, and the split between them is the point.** The CPU
reference sorts each window's exact 256-bit keys and is the arbiter; the device
inserts them into an `atomicCAS` open-addressing table at ~1e9 pairs/s and is a
*filter* whose every candidate the host re-checks exactly. They share the window
partition (which is the definition of which pairs exist, tested against brute
force) and share nothing else.

**Correctness, given that S-E has no oracle.** Nothing catalogues champion
collisions, and at $D = 2^{19}$ there is no brute force, so the fixtures are:
`champion_meet` agreeing with the *production* meet and `meet_oracle` on real
small-D point sets, where all three are valid; a planted collision at full width
(65 = 1²+8² = 4²+7², shifted so ν is ~122 bits and the key ~244) that must come
back; and `stats.pairs == C(D,2)` on every run, which is the statement that no
pair was skipped — the failure mode that has no other symptom.

## The local model (G1)

`./build/dev/apps/tower_demo [max_modulus] [max_depth] [only_p]` computes
**V_p(k,d)**: the least cap on node valuations admitting a depth-*d* tower over
ℤ/pᵏ. The default budget reaches 3⁷ = 2187 in about ten seconds.

**Why mod pᵏ.** [`drafts/gpu_implementation_plan.md`](../../drafts/gpu_implementation_plan.md)
§2.2 rejects two candidate local models — "roots stay distinct mod *p*" and
"splits mod *p* at all" — and leaves σ_p undefined. They are not two competing
answers. A solution is a ℤ_p-point of the *nondegenerate* locus, and
nondegenerate means every node value has **finite** valuation; mod *p* you can
see zero and nonzero but not finite. Bounding the valuation instead recovers
both as the two extreme caps at k = 1: cap 0 is the distinct-root model, cap k
is the multiplicity model. `tower_test.cc` asserts exactly that, against the
numbers `drafts/verify_gpu_plan.py` prints for each.

**What a stable entry means.** If V_p(k,d) < k then every E_d over ℤ has a node
value divisible by p^V — unconditionally, since a real tower reduces to a tower
mod pᵏ and reduction can only raise a valuation to the cap. Measured, this
agrees with `paper.md` Theorem 4.1 everywhere and exceeds it in several places;
at p = 3, d = 4 it gives 3⁶ where the theorem gives 3³, and the §7 exemplar
attains 3⁶ exactly, so that one is sharp.

**Cost, and why it is single-threaded.** One level costs `states × p^k × p^k/64`
word-operations, so the wall clock scales as `p^{3k}`: 2.6 s at 3⁷, 63 s at 3⁸,
1733 s at 3⁹ — 27.5× per step of *k*, against the predicted 27. That is the
binding constraint on §6.6, not the hardware —
and threads are the wrong lever for it. Measured at 3⁸, level 2 performs
4.3×10⁷ hash-table insertions to produce 32,806 distinct states, a **99.92% hit
rate**, because the `m` cyclic rotations of each state re-discover one
translation orbit `m` times. Threading that loop means writing a concurrent
intern table to accelerate work that is redundant by construction. Canonicalising
once per state instead removes it — a factor `m/64`, ≈300× at 3⁹ — and leaves a
per-state loop with a once-per-state critical section that threads cleanly
afterwards. See `drafts/local_model_plan.md` §5.1–§5.2. (Sweeps over different
primes are independent processes, so that parallelism is free at the shell.)

**Correctness.** `tower_oracle.cc` walks `paper.md` §2's tree of nested radicals
one L-tuple at a time — no sets, no interning, no rotations, no memoization — so
what it evaluates is the definition. `tower_count_brute` runs it over every
tuple in (ℤ/pᵏ)^depth and must equal the DP's count exactly, at every cap. It
also answers what the set DP structurally cannot: the exemplar's own tuple
recovers ht_p (6 at p = 3, 2 at p = 7) from the reduced residues alone, which is
what makes a null result at depth 5 distinguishable from insufficient precision.

## Layout

```
src/core/     neutral: no GMP, no Boost, no CUDA. Both halves depend on it.
  uint.h      pte_core           UInt<N> over 32-bit limbs, host + device
  window.h                       the S-E value-window partition
  collision.h                    the shared 4-index hit type
src/cpu/      per-app directories, GMP + Boost, never sees nvcc
  bigint/     pte_cpu_bigint     big-integer arithmetic + limb export
  matmul/     pte_cpu_matmul
  gaussian/   pte_cpu_gaussian   Gaussian arithmetic, Cornacchia, point sets
  enumerate/  pte_cpu_enumerate  norm generation, reference meet, primitivity
  champion/   pte_cpu_champion   S-E: the reference meet and champion point sets
  tower/      pte_cpu_tower      G1: the mod-p^k tower DP and its tree oracle.
                                 Links nothing -- no GMP, no Boost, not even
                                 pte_core
src/gpu/      per-app directories, cudart only, never sees GMP
  core/       pte_gpu_core       error checking, launch macro, host/device
                                 buffers, pinned pool, device discovery.
                                 No kernels.
  bigint/     pte_gpu_bigint
  matmul/     pte_gpu_matmul
  meet/       pte_gpu_meet       the E4 meet
  champion/   pte_gpu_champion   S-E: the hash-insert kernel
src/bindings/ pte_native._core   pybind11 module linking both
apps/         enum_demo          the E4 enumeration, CPU+GPU
              champion_demo      S-E, CPU+GPU
              tower_demo         G1, CPU only
              hybrid_demo        the original scaffold demo
python/pte_native/               the Python package and its tests
```

`src/core/` is the one place a header is compiled by *both* compilers. It exists
because S-E's key crosses the PCIe bus as raw bytes, so the two halves have to
agree on a layout by construction rather than by review; it is kept free of GMP,
Boost and the CUDA runtime so that depending on it creates no edge in either
direction. See REVIEW.md D5, which anticipated it.

Each app directory owns its library *and* its test binary, so a suite links
only what it covers. `pte_cpu` and `pte_gpu` are INTERFACE targets aggregating
their halves, so consumers can still link "the CPU half" without naming apps.

The include root is `src/`, so every header is reached by its full path —
`#include "gpu/meet/meet.h"`, `#include "cpu/gaussian/gaussian.h"`. That is not
cosmetic: `cpu/bigint/bigint.h` and `gpu/bigint/bigint.h` are different headers
with the same basename, and a bare `#include "bigint.h"` would silently resolve
by include-path order.

Dependencies point one way only: apps depend on `core`, never on each other.

Tests live next to the code they cover (`bigint.cc` ↔ `bigint_test.cc`).

### The seam: the per-app `*.h`

Each GPU app exposes one plain-C++ header — [`gpu/meet/meet.h`](src/gpu/meet/meet.h),
[`gpu/bigint/bigint.h`](src/gpu/bigint/bigint.h),
[`gpu/matmul/matmul.h`](src/gpu/matmul/matmul.h),
[`gpu/core/device.h`](src/gpu/core/device.h) — with no `<cuda_runtime.h>`, no
`__global__`, and no CUDA types in any signature. Everything outside a `.cu` is
therefore compiled by `g++`, including the demos, the bindings and the CPU tests.

### The file-extension rule

> **`.cu`/`.cuh` iff the file contains device code** — a `__device__`/`__global__`
> definition, `<<< >>>`, or a device intrinsic like `atomicAdd`. Everything else
> is `.cc`/`.h`, *including every call to the CUDA runtime API*.

`cudaMalloc`, `cudaMemcpy`, `cudaHostAlloc` and friends are a C library: they
need `<cuda_runtime.h>` and `-lcudart`, not nvcc. So `gpu/core/check.h`,
`device_buffer.h`, `host_buffer.h` and `device_append_buffer.h` are all ordinary
headers, and only three files in `src/gpu` are `.cu` at all.

The rule is load-bearing rather than cosmetic: **an unguarded `__device__`
function in a header promotes every translation unit that includes it,
transitively, to nvcc.** `PTE_CUDA_LAUNCH_KERNEL` is therefore split out of
`check.h` into [`launch.cuh`](src/gpu/core/launch.cuh) — `<<< >>>` is the one
thing in it that genuinely needs nvcc.

The escape hatch, where a type wants both a host owner and a device method, is
`#ifdef __CUDACC__` around the method:
[`DeviceAppendView::push`](src/gpu/core/device_append_buffer.h) is a `__device__`
member that host translation units simply do not see, so the header stays a
plain `.h` and `meet.cc` stays ordinary C++. This gives the class different
definitions in host and device TUs — formally an ODR violation, benign here
because the members are identical (layout and the by-value kernel argument are
unaffected) and `push` is inline `__device__`, emitted only into device code.
It is the standard CUDA idiom; Thrust and CUB are built on it. **Anything added
inside such a guard must not touch layout.**

The same hatch, widened for a whole value type, is `PTE_HOST_DEVICE` in
[`src/core/hostdev.h`](src/core/hostdev.h): empty under `g++`, `__host__
__device__` under `nvcc`. `core/uint.h` needs it because `UInt<N>` is arithmetic
the host and the device both perform on the same bytes. The obligation it
creates is checked rather than asserted —
[`gpu/bigint/uint_device_test.cu`](src/gpu/bigint/uint_device_test.cu) runs a
carry-heavy corpus through both definitions and requires identical results,
because `add` and `sub` genuinely have two bodies there (a PTX carry chain on
the device, a portable accumulate on the host). Dispatching on `__CUDA_ARCH__`
to reach a faster instruction for the *same* value is fine and is what that test
polices; dispatching to a different value is not.

Practical payoff: host-side files get `-Wpedantic`, which CUDA translation units
in this project do not.

That one constraint buys three things:

1. Each test binary links only the app it covers, so the suites are *provably*
   independent — the GPU tests write their own host oracles rather than calling
   GMP, and nothing in `src/gpu` can reach `pte_cpu`.
2. `cmake --preset cpu-only` works on a machine with no CUDA toolkit at all.
3. Adding a kernel never drags nvcc into the rest of the build.

An app with a kernel follows a three-layer shape, `meet` being the worked
example:

| Layer | File | Sees | Compiler |
|---|---|---|---|
| public API | `meet.h` | `std::span`, plain C++ only | — |
| host side | `meet.cc` | validation, staging, device buffers, launch | `g++` |
| device contract | `meet_impl.h` | `DeviceBuffer<T>`, the append view | — |
| kernel | `meet_impl.cu` | the kernel itself | `nvcc` |

Neither API takes a `(pointer, count)` pair. The public one takes `std::span`
and the device contract takes `const DeviceBuffer<T>&`; both carry their own
extent, so a caller cannot pass a length that disagrees with its data — which
had been the substance of most of the validation. **Raw device pointers appear
exactly once in the app, in the kernel's argument list**, which cannot take
anything else.

Tests go through the public header, so they cover validation and staging rather
than just the kernel — and `meet_test.cc` being a `.cc` is mechanical proof of
the seam: a test that cannot contain device code cannot reach past `meet.h`.

The contract crossing the big-integer seam is the limb layout: **a big integer
is a fixed-width array of 32-bit limbs, least-significant limb first.**
`pte::cpu::to_limbs` / `from_limbs` convert to and from `mpz_class`.

## Presets

| Preset | Build type | CUDA | Notes |
|---|---|---|---|
| `dev` | RelWithDebInfo | on | default; `-lineinfo` for profiling |
| `debug` | Debug | on | |
| `release` | Release | on | |
| `cpu-only` | RelWithDebInfo | **off** | builds and tests without any CUDA toolkit |

All CUDA presets target `sm_89` (RTX 4060, Ada). Change
`CMAKE_CUDA_ARCHITECTURES` if you move to another GPU.

## Testing

CTest labels let the suites run separately:

```sh
ctest --preset dev -L cpu   # no GPU needed
ctest --preset dev -L gpu
```

The GPU tests also `GTEST_SKIP` themselves when no device is present, so they
report as skipped rather than failed on a GPU-less machine.

Oracles are always a second, independent implementation:

- CPU big integers → `boost::multiprecision::cpp_int`
- `core::UInt<N>` → `boost::multiprecision::uint256_t` on the host, and its own
  host definition when run on the device (`uint_device_test.cu`)
- CPU matmul → `boost::numeric::ublas::prod`
- GPU kernels → host loops written inline in the test files
- the champion meet → three ways: the production meet and `meet_oracle` at small
  D, a `std::map` brute force inside the GPU test, and a planted collision at
  full width
- Python bindings → Python's own `int` and numpy's `@`

Kernel memory and race correctness:

```sh
compute-sanitizer --tool memcheck  ./build/dev/src/gpu/meet/pte_gpu_meet_test
compute-sanitizer --tool racecheck ./build/dev/src/gpu/champion/pte_gpu_champion_test
compute-sanitizer --tool initcheck ./build/dev/src/gpu/champion/pte_gpu_champion_test
```

`racecheck` matters most on `gpu/champion`, which is the only kernel here whose
correctness is an argument about concurrency rather than about arithmetic: its
hash insert relies on "of two equal keys, exactly one sees the other", which is
a claim about `atomicCAS` ordering. All three tools are clean on it — but note
what they cannot see. The false-sharing class of bug (per-thread state sharing a
cache line) is invisible to ThreadSanitizer on the CPU and its device analogues
— shared-memory bank conflicts, contention on an append counter — are invisible
to `compute-sanitizer`. Those need `ncu` sector and conflict counters, which are
[currently blocked on this machine](../../drafts/performance_ceiling.md).

## Python

```sh
pip install -e .   # scikit-build-core drives CMake
pytest
```

```python
import pte_native as pn

pn.mul(2**500, 3**300)          # CPU, via GMP; exact Python int in and out
pn.mul_gpu(2**200, 3**150)      # CPU packs limbs → CUDA kernel → CPU unpacks
pn.mul_batch_gpu(xs, ys)        # the batched form the GPU actually wants
pn.matmul(a, b)                 # 2-D float64 numpy arrays
pn.matmul_gpu(a, b)
pn.cuda_available(), pn.device_name(0)
```

Integers convert through `mpz_class` via the caster in
[`src/bindings/mpz_caster.h`](src/bindings/mpz_caster.h), so arbitrary-precision
values survive the round trip exactly. The `_gpu` functions always exist; in a
CPU-only build they raise `RuntimeError`, and `pn.HAVE_CUDA` is `False`.

Three things to know:

- **After editing C++ or CUDA, re-run `pip install -e .`** to pick the change up.
  scikit-build-core can rebuild automatically on import instead, but only with
  build isolation disabled — with isolation on, the ephemeral
  `/tmp/pip-build-env-*/…/cmake` path gets baked into `build/pip/build.ninja`,
  the install succeeds, and every later import then fails. So if you want the
  fast loop, opt in explicitly and pass both flags:

  ```sh
  pip install -e . --no-build-isolation --config-settings=editable.rebuild=true
  ```

  (That needs `scikit-build-core`, `pybind11`, `cmake` and `ninja` in the venv,
  since nothing is fetched for you.)

- **`build/pip` is a persistent CMake cache.** Passing a different
  `--config-settings=cmake.define.X=Y` does *not* invalidate it. To change a
  CMake option for the pip build, `rm -rf build/pip` first. For example, to get a
  CPU-only module:

  ```sh
  rm -rf build/pip
  pip install -e . --config-settings=cmake.define.PTE_ENABLE_CUDA=OFF
  ```

- Tests are named `*_test.py`, not `test_*.py`; `pyproject.toml` configures
  pytest accordingly.

## Dependency discovery, and why it looks like this

- **Boost** is found in CONFIG mode (`find_package(Boost CONFIG)` →
  `Boost::headers`). Module-mode `FindBoost` is deprecated and this builds under
  CMake 4. Only header-only Boost is used, because this machine has no compiled
  Boost libraries installed (no `libboost_unit_test_framework`, hence GoogleTest
  rather than Boost.Test).
- **GMP** has no CMake find module, but does ship `gmp.pc` / `gmpxx.pc`, so it
  comes in via `pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp gmpxx)`. No
  hand-written `FindGMP.cmake`.
- **GoogleTest** is pinned to `v1.17.0` via `FetchContent` (needs network on
  first configure). 1.14 declares `cmake_minimum_required(VERSION 3.13)`, which
  CMake 4 rejects outright.

## Sample output

```
$ ./build/dev/apps/enum_demo
CUDA devices: 1
Using device 0: NVIDIA GeForce RTX 4060 Laptop GPU
Device warm-up (context creation): 53.08 ms

== E_4 enumeration, admissible L1 <= 2000000 ==
  admissible norms        23576
  points                 108909
  pairs                  219812

  CPU norm generation     1.29 ms
  CPU point sets         12.12 ms
  CPU pack to CSR         1.25 ms
  GPU meet                0.94 ms  (includes H2D+D2H copies)
  CPU verify hits         0.00 ms
  CPU reference meet      1.27 ms

  raw GPU hits                 3
    rejected, key mismatch in wide arithmetic 0
    rejected, imprimitive (scaled duplicate)  1
  primitive E_4 solutions: GPU 2, CPU 2
  CPU and GPU solution sets are identical -> OK
    E_4 at L1 = 67405: nu = 4037 27643 35653 44933
    E_4 at L1 = 600445: nu = 151093 332533 515573 594613

  transferred 0.83 MB to do ~219812 pair comparisons
  arithmetic intensity ~0.25 ops/byte -- transfer bound;
  generating norms on the device is what fixes this.

All checks passed.
```

Both solutions check out independently: $L_1 = 67405 = 5\cdot13\cdot17\cdot61$
has root pairs $(11,367), (77,359), (101,353), (131,343)$, and the recovered
tower satisfies the degree-16 identity
$\prod(x^2-r_k^2) = (((x^2-L_1)^2-L_2)^2-L_3)^2-L_4$ exactly. Note the third raw
hit being rejected as a scaled duplicate — that is the primitivity rule doing
real work, not a theoretical concern.

The champion meet, at the sized champion (this is the full S-E search over one
norm — 1.374×10¹¹ pairs; the default `champion_demo` with no argument is the
14-prime version and takes three seconds):

```
$ ./build/dev/apps/champion_demo 20
== S-E champion meet, 20 split primes ==
  L1                     2386442377984630986454048784862092465
  L1 bits                121
  points D               524288
  pairs C(D,2)           1.374e+11
  widest nu / key        121 / 243 bits
  point-set construction 1.009 s   (GMP; the int32 path would overflow above 20-bit coordinates)

  CPU reference (sort each window, exact 256-bit keys), D = 8192
    windows 42, pairs 33550336, hits 0
    2.306 s   1.45e+07 pair/s

  hash table 2^26 slots (512 MB), window cap 3.36e+07 pairs

  GPU (atomicCAS hash insert per window, host confirms)
    windows 5794, pairs 137438691328, candidates 4836, false positives 4836, hits 0
    278.060 s total   4.94e+08 pair/s   (11588 H2D copies)
      device 138.776 s  ->  9.9e+08 pair/s in-kernel
      host window planning 139.252 s  ->  50% of the run, device idle throughout

  cross-check at D = 8192: CPU 2.306 s, GPU 0.043 s (54x)
  CPU and GPU hit sets are identical -> OK  (0 hits)

  planted collision (nu ~ 2^126, key ~ 2^253): CPU found, GPU found, sets identical

All checks passed.
```

Read the middle rows rather than the total. `pairs 137438691328` is
$\binom{2^{19}}{2}$ **exactly** — every pair enumerated once, which is the
property value-windowing exists to have and the one whose failure is invisible.
All 4,836 candidates were fingerprint false positives and all were rejected by
the exact re-check, at the rate the model predicts (~0.8 per window against ~0.5
expected). This particular champion has no E₄; that is a fact about one norm and
excludes nothing, since S-E is a discovery search.

The last two rows are the ones to read. The device half, 138.8 s at 9.9e8
pair/s, **is** bench P4's measured `atomicCAS` ceiling of 1.0e9 insert/s — the
cost model projected 138 s from a primitive benchmark and the search took 138.8.
The kernel has nothing left to give at this table size. The other half is the
host planning an O(D) window sweep 5,794 times while the device idles, and that
is where the next work is: see plan §3.4 for the ordering (overlap the planner,
then shrink the table into L2 for a measured 2.9×, then partition the hash).

The original scaffold demo:

```
$ ./build/dev/apps/hybrid_demo
CUDA devices: 1
Using device 0: NVIDIA GeForce RTX 4060 Laptop GPU
Device warm-up (context creation): 58.55 ms

== Big integers: 20000 pairs of 256-bit operands ==
  CPU (GMP mul)          1.27 ms
  CPU pack to limbs      1.84 ms
  GPU (schoolbook)       1.79 ms  (includes H2D+D2H copies)
  CPU verify vs GMP      1.14 ms
  mismatches: 0 / 20000 -> OK

== Matrix multiply: (512x512) * (512x512), double ==
  CPU (i-k-j loops)     34.07 ms  (7.88 GFLOP/s)
  GPU (16x16 tiled)      4.72 ms  (56.89 GFLOP/s, incl. copies)
  max |cpu - gpu|: 2.132e-14 -> OK

All checks passed.
```

Note what this does *not* show: the naive schoolbook kernel is only about level
with GMP at 256 bits, because one thread per multiplication leaves the GPU
mostly idle and the host↔device copies dominate. The matmul, which has real
arithmetic intensity, is ~7× faster than the CPU. Creating the CUDA context
costs ~60 ms and happens on the first launch, which is why the demo warms up
before timing anything.

## Caveats

- The GPU big-integer kernels are **unsigned** and capped at
  `pte::gpu::kMaxLimbs` = 32 limbs (1024 bits), since the accumulator is a
  per-thread local array. `mul_gpu` raises on negative or oversized operands.
- `double` on a 4060 runs at 1/64 of `float` rate. Fine at these sizes, and it
  keeps the CPU/GPU comparison exact enough to assert on.
- CUDA 12.9 binaries run on this driver (550 / CUDA 12.4) thanks to 12.x minor
  version compatibility. If that ever breaks, the fallback is `g++-12` with
  `/usr/bin/nvcc`; both are installed.
