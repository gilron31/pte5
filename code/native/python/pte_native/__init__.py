"""Python bindings for the pte5 C++/CUDA sandbox.

Two paths are exposed from one import:

* CPU only, backed by GMP and Boost -- :func:`add`, :func:`mul`, :func:`matmul`.
* CPU + GPU, where the host packs the data and a CUDA kernel does the
  arithmetic -- :func:`mul_gpu`, :func:`mul_batch_gpu`, :func:`matmul_gpu`.

The integer functions take and return native Python ints of arbitrary size.
The matrix functions take and return 2-D float64 numpy arrays.

    >>> import pte_native as pn
    >>> pn.mul(2**200, 3**150) == 2**200 * 3**150
    True
    >>> pn.cuda_available() and pn.mul_gpu(2**200, 3**150) == 2**200 * 3**150
    True

``HAVE_CUDA`` says whether this build has the CUDA half compiled in at all;
:func:`cuda_available` says whether a usable device is present right now. The
``_gpu`` functions always exist, and raise ``RuntimeError`` when ``HAVE_CUDA``
is False.
"""

from ._core import (
    HAVE_CUDA,
    MAX_GPU_LIMBS,
    add,
    cuda_available,
    device_count,
    device_name,
    matmul,
    matmul_gpu,
    mul,
    mul_batch_gpu,
    mul_gpu,
)

__all__ = [
    "HAVE_CUDA",
    "MAX_GPU_LIMBS",
    "add",
    "cuda_available",
    "device_count",
    "device_name",
    "matmul",
    "matmul_gpu",
    "mul",
    "mul_batch_gpu",
    "mul_gpu",
]
