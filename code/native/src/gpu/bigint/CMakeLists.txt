add_library(pte_gpu_bigint STATIC bigint_kernel.cu)
target_link_libraries(pte_gpu_bigint PUBLIC pte_gpu_core PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  # Two suites in one binary:
  #   bigint_kernel_test  the schoolbook multiply kernel (arbitrary width)
  #   uint_device_test    core/uint.h's device path against its host path --
  #                       the check that makes PTE_HOST_DEVICE safe, and the
  #                       reason this app is the one the S-E plan promotes
  add_executable(pte_gpu_bigint_test bigint_kernel_test.cu uint_device_test.cu)
  target_link_libraries(pte_gpu_bigint_test
    PRIVATE pte_gpu_bigint pte_core GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_bigint_test COMMAND pte_gpu_bigint_test)
  set_tests_properties(pte_gpu_bigint_test PROPERTIES LABELS gpu)
endif()
