add_library(pte_gpu_matmul STATIC matmul_kernel.cu)
target_link_libraries(pte_gpu_matmul PUBLIC pte_gpu_core PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_gpu_matmul_test matmul_kernel_test.cu)
  target_link_libraries(pte_gpu_matmul_test
    PRIVATE pte_gpu_matmul GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_matmul_test COMMAND pte_gpu_matmul_test)
  set_tests_properties(pte_gpu_matmul_test PROPERTIES LABELS gpu)
endif()
