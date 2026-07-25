// The `pte_native._core` extension module.
//
// Exposes both paths asked of this sandbox: pure-CPU entry points backed by GMP
// and Boost, and CPU+GPU entry points that pack operands on the host, run the
// CUDA kernel, and unpack the result. The `_gpu` functions exist unconditionally
// so calling code can rely on them being present; they raise RuntimeError when
// the module was built without CUDA.
#include "mpz_caster.h" // must precede any binding that mentions mpz_class

#include "bigint.h"
#include "matmul.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // list <-> std::vector for the batched entry points

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef PTE_HAVE_CUDA
#include "gpu.h"
#endif

namespace py = pybind11;

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

// Shape checking shared by the CPU and GPU matmul bindings.
struct MatmulShape {
  int M, N, K;
};

MatmulShape CheckMatmulShapes(const Array &A, const Array &B) {
  if (A.ndim() != 2 || B.ndim() != 2) {
    throw std::invalid_argument("matmul: both operands must be 2-D, got " +
                               std::to_string(A.ndim()) + "-D and " +
                               std::to_string(B.ndim()) + "-D");
  }
  const auto M = A.shape(0), K = A.shape(1);
  const auto K2 = B.shape(0), N = B.shape(1);
  if (K != K2) {
    throw std::invalid_argument(
        "matmul: shape mismatch, (" + std::to_string(M) + "," +
        std::to_string(K) + ") @ (" + std::to_string(K2) + "," +
        std::to_string(N) + ")");
  }
  return {static_cast<int>(M), static_cast<int>(N), static_cast<int>(K)};
}

Array MatmulCpu(const Array &A, const Array &B) {
  const auto s = CheckMatmulShapes(A, B);
  Array C({s.M, s.N});
  const double *a = A.data();
  const double *b = B.data();
  double *c = C.mutable_data();
  {
    py::gil_scoped_release unlock;
    pte::cpu::matmul(a, b, c, s.M, s.N, s.K);
  }
  return C;
}

#ifdef PTE_HAVE_CUDA

Array MatmulGpu(const Array &A, const Array &B) {
  const auto s = CheckMatmulShapes(A, B);
  Array C({s.M, s.N});
  const double *a = A.data();
  const double *b = B.data();
  double *c = C.mutable_data();
  {
    py::gil_scoped_release unlock;
    pte::gpu::matmul(a, b, c, s.M, s.N, s.K);
  }
  return C;
}

// The hybrid path for one pair: GMP -> limbs -> CUDA kernel -> limbs -> GMP.
mpz_class MulGpu(const mpz_class &a, const mpz_class &b) {
  if (a < 0 || b < 0) {
    throw std::invalid_argument(
        "mul_gpu: the GPU kernels only handle non-negative operands");
  }
  const std::size_t needed =
      std::max<std::size_t>({pte::cpu::limb_count(a), pte::cpu::limb_count(b), 1});
  if (needed > static_cast<std::size_t>(pte::gpu::kMaxLimbs)) {
    throw std::invalid_argument(
        "mul_gpu: operands need " + std::to_string(needed) +
        " limbs, kernel maximum is " + std::to_string(pte::gpu::kMaxLimbs) +
        " (" + std::to_string(32 * pte::gpu::kMaxLimbs) + " bits)");
  }
  const int nlimbs = static_cast<int>(needed);
  const auto limbs_a = pte::cpu::to_limbs(a, nlimbs);
  const auto limbs_b = pte::cpu::to_limbs(b, nlimbs);
  std::vector<std::uint32_t> out(static_cast<std::size_t>(2) * nlimbs);
  {
    py::gil_scoped_release unlock;
    pte::gpu::mul_batch(limbs_a.data(), limbs_b.data(), out.data(), nlimbs, 1);
  }
  return pte::cpu::from_limbs(out.data(), out.size());
}

// Batched variant: the whole point of the GPU path is doing many at once.
std::vector<mpz_class> MulBatchGpu(const std::vector<mpz_class> &a,
                                   const std::vector<mpz_class> &b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("mul_batch_gpu: operand lists must be the same "
                                "length, got " +
                                std::to_string(a.size()) + " and " +
                                std::to_string(b.size()));
  }
  if (a.empty())
    return {};

  // One common width for the whole batch, driven by the widest operand.
  std::size_t needed = 1;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] < 0 || b[i] < 0) {
      throw std::invalid_argument("mul_batch_gpu: the GPU kernels only handle "
                                  "non-negative operands (index " +
                                  std::to_string(i) + ")");
    }
    needed = std::max({needed, pte::cpu::limb_count(a[i]),
                       pte::cpu::limb_count(b[i])});
  }
  if (needed > static_cast<std::size_t>(pte::gpu::kMaxLimbs)) {
    throw std::invalid_argument(
        "mul_batch_gpu: operands need " + std::to_string(needed) +
        " limbs, kernel maximum is " + std::to_string(pte::gpu::kMaxLimbs));
  }

  const int nlimbs = static_cast<int>(needed);
  const int batch = static_cast<int>(a.size());
  std::vector<std::uint32_t> flat_a(static_cast<std::size_t>(batch) * nlimbs);
  std::vector<std::uint32_t> flat_b(flat_a.size());
  for (int t = 0; t < batch; ++t) {
    const std::size_t off = static_cast<std::size_t>(t) * nlimbs;
    pte::cpu::to_limbs(a[t], flat_a.data() + off, nlimbs);
    pte::cpu::to_limbs(b[t], flat_b.data() + off, nlimbs);
  }

  std::vector<std::uint32_t> flat_out(flat_a.size() * 2);
  {
    py::gil_scoped_release unlock;
    pte::gpu::mul_batch(flat_a.data(), flat_b.data(), flat_out.data(), nlimbs,
                        batch);
  }

  std::vector<mpz_class> out(batch);
  for (int t = 0; t < batch; ++t) {
    out[t] = pte::cpu::from_limbs(
        flat_out.data() + static_cast<std::size_t>(t) * 2 * nlimbs, 2 * nlimbs);
  }
  return out;
}

#else // !PTE_HAVE_CUDA

[[noreturn]] void NoCuda(const char *what) {
  throw std::runtime_error(std::string(what) +
                           ": this build of pte_native was configured without "
                           "CUDA (PTE_ENABLE_CUDA=OFF)");
}

#endif // PTE_HAVE_CUDA

} // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "CPU (GMP/Boost) and CPU+GPU (CUDA) kernels for the pte5 sandbox.";

#ifdef PTE_HAVE_CUDA
  m.attr("HAVE_CUDA") = true;
  m.attr("MAX_GPU_LIMBS") = pte::gpu::kMaxLimbs;
#else
  m.attr("HAVE_CUDA") = false;
  m.attr("MAX_GPU_LIMBS") = 0;
#endif

  // --- CPU: big integers, backed by GMP -------------------------------------
  m.def("add", &pte::cpu::add, py::arg("a"), py::arg("b"),
        "Add two integers on the CPU using GMP.");
  m.def("mul", &pte::cpu::mul, py::arg("a"), py::arg("b"),
        "Multiply two integers on the CPU using GMP.");

  // --- CPU: matrix multiply -------------------------------------------------
  m.def("matmul", &MatmulCpu, py::arg("a"), py::arg("b"),
        "Multiply two 2-D float64 arrays on the CPU.");

  // --- CPU+GPU --------------------------------------------------------------
#ifdef PTE_HAVE_CUDA
  m.def("cuda_available", []() { return pte::gpu::device_count() > 0; },
        "True if a usable CUDA device is present.");
  m.def("device_count", &pte::gpu::device_count, "Number of CUDA devices.");
  m.def("device_name", &pte::gpu::device_name, py::arg("index") = 0,
        "Name of the given CUDA device.");
  m.def("mul_gpu", &MulGpu, py::arg("a"), py::arg("b"),
        "Multiply two non-negative integers with the CUDA schoolbook kernel.");
  m.def("mul_batch_gpu", &MulBatchGpu, py::arg("a"), py::arg("b"),
        "Multiply two equal-length lists of non-negative integers on the GPU.");
  m.def("matmul_gpu", &MatmulGpu, py::arg("a"), py::arg("b"),
        "Multiply two 2-D float64 arrays with the tiled CUDA kernel.");
#else
  m.def("cuda_available", []() { return false; },
        "True if a usable CUDA device is present.");
  m.def("device_count", []() { return 0; }, "Number of CUDA devices.");
  m.def("device_name", [](int) -> std::string { NoCuda("device_name"); },
        py::arg("index") = 0, "Name of the given CUDA device.");
  m.def("mul_gpu",
        [](const mpz_class &, const mpz_class &) -> mpz_class {
          NoCuda("mul_gpu");
        },
        py::arg("a"), py::arg("b"), "Unavailable in this build.");
  m.def("mul_batch_gpu",
        [](const std::vector<mpz_class> &,
           const std::vector<mpz_class> &) -> std::vector<mpz_class> {
          NoCuda("mul_batch_gpu");
        },
        py::arg("a"), py::arg("b"), "Unavailable in this build.");
  m.def("matmul_gpu",
        [](const Array &, const Array &) -> Array { NoCuda("matmul_gpu"); },
        py::arg("a"), py::arg("b"), "Unavailable in this build.");
#endif
}
