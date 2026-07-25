"""Tests for the Python bindings.

Python's own arbitrary-precision int is the oracle for the big-integer paths and
numpy's ``@`` for the matrix paths, so these tests check the bindings (and the
int <-> mpz_class conversion) rather than re-checking the C++ that already has
its own suites.
"""

import numpy as np
import pytest

import pte_native as pn

requires_cuda = pytest.mark.skipif(
    not pn.cuda_available(), reason="no usable CUDA device"
)

# Values chosen to straddle 32-bit limb boundaries and the int/long divide.
INTS = [
    0,
    1,
    2,
    2**31 - 1,
    2**32 - 1,
    2**32,
    2**64 - 1,
    2**64,
    2**127 - 1,
    2**128,
    10**30,
    123456789012345678901234567890,
]


class TestCpuBigInt:
    @pytest.mark.parametrize("a", INTS)
    @pytest.mark.parametrize("b", INTS)
    def test_add_matches_python(self, a, b):
        assert pn.add(a, b) == a + b

    @pytest.mark.parametrize("a", INTS)
    @pytest.mark.parametrize("b", INTS)
    def test_mul_matches_python(self, a, b):
        assert pn.mul(a, b) == a * b

    def test_results_are_exact_python_ints(self):
        got = pn.mul(2**500, 2**500)
        assert isinstance(got, int)
        assert got == 2**1000  # no float rounding anywhere in the round trip

    def test_handles_negative_operands(self):
        assert pn.add(-5, 3) == -2
        assert pn.mul(-7, 6) == -42
        assert pn.mul(-(2**100), -(2**100)) == 2**200

    def test_huge_operands(self):
        a, b = 7**2000, 11**1500
        assert pn.mul(a, b) == a * b

    def test_accepts_numpy_integers(self):
        # The caster goes through __index__, so numpy scalars work.
        assert pn.mul(np.int64(6), np.int64(7)) == 42

    def test_rejects_non_integers(self):
        with pytest.raises(TypeError):
            pn.mul(1.5, 2)
        with pytest.raises(TypeError):
            pn.mul("3", 2)


class TestCpuMatmul:
    @pytest.mark.parametrize(
        "m,k,n", [(1, 1, 1), (2, 3, 4), (16, 16, 16), (37, 41, 53), (1, 64, 1)]
    )
    def test_matches_numpy(self, m, k, n):
        rng = np.random.default_rng(m * 1000 + k * 10 + n)
        a = rng.standard_normal((m, k))
        b = rng.standard_normal((k, n))
        got = pn.matmul(a, b)
        assert got.shape == (m, n)
        assert got.dtype == np.float64
        np.testing.assert_allclose(got, a @ b, rtol=1e-12, atol=1e-12)

    def test_accepts_non_contiguous_and_wrong_dtype(self):
        # forcecast in the binding means a transposed view or an int array is
        # converted rather than rejected.
        a = np.arange(6, dtype=np.int32).reshape(2, 3)
        b = np.asfortranarray(np.arange(12, dtype=np.float32).reshape(3, 4))
        np.testing.assert_allclose(pn.matmul(a, b), a @ b)

    def test_rejects_shape_mismatch(self):
        a = np.ones((2, 3))
        b = np.ones((4, 5))
        with pytest.raises(ValueError):
            pn.matmul(a, b)

    def test_rejects_non_2d(self):
        with pytest.raises(ValueError):
            pn.matmul(np.ones(3), np.ones((3, 3)))


class TestGpuAvailability:
    def test_flags_are_consistent(self):
        assert isinstance(pn.HAVE_CUDA, bool)
        if not pn.HAVE_CUDA:
            assert pn.cuda_available() is False
            assert pn.device_count() == 0
        # A device can only be usable if the CUDA half was compiled in.
        assert not pn.cuda_available() or pn.HAVE_CUDA

    @requires_cuda
    def test_device_is_nameable(self):
        assert pn.device_count() >= 1
        assert pn.device_name(0)


class TestGpuBigInt:
    @requires_cuda
    @pytest.mark.parametrize("a", INTS)
    @pytest.mark.parametrize("b", INTS)
    def test_mul_gpu_matches_python(self, a, b):
        assert pn.mul_gpu(a, b) == a * b

    @requires_cuda
    def test_mul_gpu_agrees_with_cpu_path(self):
        a, b = 2**255 - 12345, 3**160
        assert pn.mul_gpu(a, b) == pn.mul(a, b)

    @requires_cuda
    def test_mul_batch_gpu_matches_python(self):
        rng = np.random.default_rng(7)
        a = [int.from_bytes(rng.bytes(32), "little") for _ in range(500)]
        b = [int.from_bytes(rng.bytes(32), "little") for _ in range(500)]
        got = pn.mul_batch_gpu(a, b)
        assert got == [x * y for x, y in zip(a, b)]

    @requires_cuda
    def test_mul_batch_gpu_empty(self):
        assert pn.mul_batch_gpu([], []) == []

    @requires_cuda
    def test_mul_batch_gpu_rejects_ragged_lists(self):
        with pytest.raises(ValueError):
            pn.mul_batch_gpu([1, 2], [1])

    @requires_cuda
    def test_rejects_negative_operands(self):
        # The kernels are unsigned; the binding must refuse rather than wrap.
        with pytest.raises(ValueError):
            pn.mul_gpu(-1, 2)

    @requires_cuda
    def test_rejects_operands_wider_than_the_kernel(self):
        too_big = 2 ** (32 * pn.MAX_GPU_LIMBS)
        with pytest.raises(ValueError):
            pn.mul_gpu(too_big, 1)


class TestGpuMatmul:
    @requires_cuda
    @pytest.mark.parametrize(
        "m,k,n", [(1, 1, 1), (16, 16, 16), (37, 41, 53), (100, 1, 100), (15, 15, 15)]
    )
    def test_matches_numpy(self, m, k, n):
        rng = np.random.default_rng(m * 1000 + k * 10 + n)
        a = rng.standard_normal((m, k))
        b = rng.standard_normal((k, n))
        got = pn.matmul_gpu(a, b)
        assert got.shape == (m, n)
        np.testing.assert_allclose(got, a @ b, rtol=1e-12, atol=1e-12)

    @requires_cuda
    def test_agrees_with_cpu_path(self):
        rng = np.random.default_rng(11)
        a = rng.standard_normal((64, 48))
        b = rng.standard_normal((48, 32))
        np.testing.assert_allclose(
            pn.matmul_gpu(a, b), pn.matmul(a, b), rtol=1e-12, atol=1e-12
        )

    @requires_cuda
    def test_rejects_shape_mismatch(self):
        with pytest.raises(ValueError):
            pn.matmul_gpu(np.ones((2, 3)), np.ones((4, 5)))


@pytest.mark.skipif(pn.HAVE_CUDA, reason="only meaningful in a CPU-only build")
class TestCpuOnlyBuild:
    def test_gpu_entry_points_exist_but_raise(self):
        with pytest.raises(RuntimeError):
            pn.mul_gpu(2, 3)
        with pytest.raises(RuntimeError):
            pn.matmul_gpu(np.ones((2, 2)), np.ones((2, 2)))
