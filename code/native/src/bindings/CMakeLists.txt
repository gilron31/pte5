# The pybind11 extension module. Links both halves, so it exposes the CPU-only
# and the CPU+GPU paths from one import.
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

# Under `pip install`, scikit-build-core puts pybind11 on CMAKE_PREFIX_PATH. For
# a plain `cmake --preset dev` build, ask the interpreter where its copy lives.
if(NOT DEFINED pybind11_DIR AND NOT DEFINED pybind11_ROOT)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE _pybind11_cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pybind11_probe
    ERROR_QUIET
  )
  if(_pybind11_probe EQUAL 0)
    set(pybind11_DIR "${_pybind11_cmakedir}")
  endif()
endif()
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE module.cc)
target_link_libraries(_core PRIVATE pte_cpu pte_warnings)
if(PTE_HAVE_CUDA)
  target_link_libraries(_core PRIVATE pte_gpu)
  target_compile_definitions(_core PRIVATE PTE_HAVE_CUDA)
endif()

# Where `pip install` puts the built module. `python/pte_native/__init__.py`
# imports it as `pte_native._core`.
#
# Note that the source tree's python/pte_native/ has no _core module in it, so
# setting PYTHONPATH to point at a CMake build directory does not work: the
# source package shadows it. Import the module via `pip install -e .` instead,
# which is what the README documents and what the pytest suite expects.
install(TARGETS _core DESTINATION pte_native)
