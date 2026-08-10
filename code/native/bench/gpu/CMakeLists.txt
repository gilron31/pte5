# GPU primitive ceilings. Guarded by PTE_HAVE_CUDA from the parent.
#
# Links NVML (CUDA::nvml) because the clock has to be read at run time: this
# device idles at 210 MHz of 3105, so a result without its clock beside it is
# not a result. See gpu_harness.cuh.
add_executable(bench_gpu_primitives bench_gpu_primitives.cu)
target_include_directories(bench_gpu_primitives PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(bench_gpu_primitives
  PRIVATE CUDA::cudart CUDA::nvml pte_warnings
)
# CUB's DeviceRadixSort instantiates deep template chains; relaxed constexpr is
# already on via pte_warnings, and separable compilation is not needed.
set_target_properties(bench_gpu_primitives PROPERTIES CUDA_ARCHITECTURES 89)

if(PTE_BUILD_TESTS)
  # Smoke only: P1 alone, which allocates almost nothing and takes ~1 s. The
  # sort and hash cases want gigabytes and have no place in ctest.
  add_test(NAME bench_gpu_primitives_smoke COMMAND bench_gpu_primitives p1)
  set_tests_properties(bench_gpu_primitives_smoke PROPERTIES LABELS "bench;gpu")
endif()
