# The CUDA half. Depends on the CUDA runtime and nothing else -- in particular
# not on GMP or pte_cpu, which is what lets the GPU tests stand alone.
add_library(pte_gpu STATIC
  bigint_kernel.cu
  matmul_kernel.cu
)
target_include_directories(pte_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pte_gpu PUBLIC CUDA::cudart PRIVATE pte_warnings)
# Separable compilation is deliberately left off: no __device__ function is
# called across translation units here, so whole-program compilation is both
# faster and avoids the device-link step (which emits spurious "skipping
# incompatible libdl.a" warnings from nvlink).

if(PTE_BUILD_TESTS)
  add_executable(pte_gpu_test
    bigint_kernel_test.cu
    matmul_kernel_test.cu
  )
  target_link_libraries(pte_gpu_test
    PRIVATE pte_gpu GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_test COMMAND pte_gpu_test)
  # Labelled `gpu` so `ctest -L cpu` skips it entirely on a GPU-less machine.
  # The tests themselves also GTEST_SKIP when no device is present.
  set_tests_properties(pte_gpu_test PROPERTIES LABELS gpu)
endif()
