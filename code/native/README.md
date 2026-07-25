# `code/native` — C++/CUDA sandbox

A working scaffold for the compiled side of the E₅ work. It exists to prove out
five capabilities end to end, on two deliberately simple use cases:

| Capability | Where |
|---|---|
| C++ library against Boost and GMP | [`src/cpu/`](src/cpu/) → `pte_cpu` |
| CUDA kernels running on the local GPU | [`src/gpu/`](src/gpu/) → `pte_gpu` |
| CPU processing that offloads to the GPU | [`apps/hybrid_demo.cc`](apps/hybrid_demo.cc) |
| CPU and GPU unit tested **independently** | `*_test.cc` / `*_test.cu`, CTest labels `cpu` / `gpu` |
| Python bindings to both paths | [`src/bindings/`](src/bindings/) → `pte_native` |

Use case 1 is **big-integer add/multiply**: GMP on the CPU, naive schoolbook on
the GPU. Use case 2 is **matrix multiplication**: plain loops on the CPU, a
16×16 tiled shared-memory kernel on the GPU.

## Quick start

```sh
cd code/native

cmake --preset dev            # configure (fetches GoogleTest on first run)
cmake --build --preset dev    # build everything

ctest --preset dev -L cpu     # CPU unit tests
ctest --preset dev -L gpu     # GPU unit tests
ctest --preset dev            # both

./build/dev/apps/hybrid_demo  # the CPU→GPU→CPU demo

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

## Layout

```
src/cpu/     pte_cpu       C++ only: GMP + Boost. Never sees nvcc.
src/gpu/     pte_gpu       CUDA only: cudart. Never sees GMP.
src/bindings/ pte_native._core   pybind11 module linking both.
apps/        hybrid_demo   The two halves working together.
python/pte_native/         The Python package and its tests.
```

Tests live next to the code they cover (`bigint.cc` ↔ `bigint_test.cc`), and
each `src/*/CMakeLists.txt` owns both its library and its test binary.

### The seam: `src/gpu/gpu.h`

[`src/gpu/gpu.h`](src/gpu/gpu.h) declares only host-callable free functions in
plain C++ — no `<cuda_runtime.h>`, no `__global__`, no CUDA types in any
signature. Everything outside `src/gpu/*.cu` is therefore compiled by `g++`,
including the demo, the bindings and the CPU tests.

That one constraint buys three things:

1. The GPU test binary links `pte_gpu` alone and the CPU test binary links
   `pte_cpu` alone, so the two suites are *provably* independent — the GPU tests
   write their own host oracles rather than calling GMP.
2. `cmake --preset cpu-only` works on a machine with no CUDA toolkit at all.
3. Adding a kernel never drags nvcc into the rest of the build.

The contract crossing the seam is the limb layout: **a big integer is a
fixed-width array of 32-bit limbs, least-significant limb first.**
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
