# `code/native` — the compiled side of the E₅ work

Two things live here:

1. **The E₄ enumeration** — [`src/cpu/gaussian/`](src/cpu/gaussian/),
   [`src/cpu/enumerate/`](src/cpu/enumerate/), [`src/gpu/meet/`](src/gpu/meet/)
   and [`apps/enum_demo.cc`](apps/enum_demo.cc). This is search **S-A** of
   [`drafts/gpu_algorithms.md`](../../drafts/gpu_algorithms.md), in its first
   working form — see the section below.
2. **The original scaffold** — big-integer add/multiply (GMP vs. a naive
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

./build/dev/apps/enum_demo    # the E₄ enumeration demo
./build/dev/apps/hybrid_demo  # the original scaffold demo

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
4. Champion-regime kernel, the $\mathbb{F}_p$ tower sieve, the top-collision and
   bottom-lift joins, and the A4′ product-form investigation.

## Layout

```
src/cpu/      per-app directories, GMP + Boost, never sees nvcc
  bigint/     pte_cpu_bigint     big-integer arithmetic + limb export
  matmul/     pte_cpu_matmul
  gaussian/   pte_cpu_gaussian   Gaussian arithmetic, Cornacchia, point sets
  enumerate/  pte_cpu_enumerate  norm generation, reference meet, primitivity
src/gpu/      per-app directories, cudart only, never sees GMP
  core/       pte_gpu_core       error checking, launch macro, host/device
                                 buffers, device discovery. No kernels.
  bigint/     pte_gpu_bigint
  matmul/     pte_gpu_matmul
  meet/       pte_gpu_meet       the E4 meet
src/bindings/ pte_native._core   pybind11 module linking both
apps/         enum_demo          the E4 enumeration, CPU+GPU
              hybrid_demo        the original scaffold demo
python/pte_native/               the Python package and its tests
```

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

Anything that must name a CUDA type lives in a `.cuh` included only from `.cu`
files: `gpu/core/check.cuh`, `device_buffer.cuh`, `host_buffer.cuh`,
`device_append_buffer.cuh`, and each app's `*_impl.cuh`.

That one constraint buys three things:

1. Each test binary links only the app it covers, so the suites are *provably*
   independent — the GPU tests write their own host oracles rather than calling
   GMP, and nothing in `src/gpu` can reach `pte_cpu`.
2. `cmake --preset cpu-only` works on a machine with no CUDA toolkit at all.
3. Adding a kernel never drags nvcc into the rest of the build.

An app with a kernel follows a three-layer shape, `meet` being the worked
example:

| Layer | File | Sees |
|---|---|---|
| public API | `meet.h` | plain C++ only |
| host side | `meet.cu` | validation, staging, device buffers, launch |
| device contract | `meet_impl.cuh` | device pointers |
| kernel | `meet_impl.cu` | the kernel itself |

Tests go through the public header, so they cover validation and staging rather
than just the kernel.

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
- CPU matmul → `boost::numeric::ublas::prod`
- GPU kernels → host loops written inline in the test files
- Python bindings → Python's own `int` and numpy's `@`

Kernel memory and race correctness:

```sh
compute-sanitizer --tool memcheck  ./build/dev/src/gpu/pte_gpu_test
compute-sanitizer --tool racecheck ./build/dev/src/gpu/pte_gpu_test --gtest_filter='GpuMatmul.*'
```

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
