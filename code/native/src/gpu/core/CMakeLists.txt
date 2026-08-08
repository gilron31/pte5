# Shared GPU infrastructure: error checking, the launch macro, the owning
# host/device buffers, and device discovery. No kernels of its own.
add_library(pte_gpu_core STATIC device.cu)
# Include root is src/, so every header is reached as "gpu/core/...", which
# keeps gpu/bigint/bigint.h distinct from cpu/bigint/bigint.h.
target_include_directories(pte_gpu_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_gpu_core PUBLIC CUDA::cudart PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_gpu_core_test
    device_append_buffer_test.cu
    host_buffer_test.cu
  )
  target_link_libraries(pte_gpu_core_test
    PRIVATE pte_gpu_core GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_core_test COMMAND pte_gpu_core_test)
  set_tests_properties(pte_gpu_core_test PROPERTIES LABELS gpu)
endif()
