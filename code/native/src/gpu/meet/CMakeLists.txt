# The E_4 meet, in the three-layer shape:
#   meet.h        plain C++ API (no CUDA types), so g++ TUs can call it
#   meet.cu       host side: validation, staging, device buffers, launch
#   meet_impl.cuh device-side contract over device pointers
#   meet_impl.cu  the kernel
add_library(pte_gpu_meet STATIC meet.cu meet_impl.cu)
target_link_libraries(pte_gpu_meet PUBLIC pte_gpu_core PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  # Goes through meet.h, so the test covers validation and staging as well as
  # the kernel, and links pte_gpu_meet alone -- never pte_cpu.
  add_executable(pte_gpu_meet_test meet_test.cu)
  target_link_libraries(pte_gpu_meet_test
    PRIVATE pte_gpu_meet GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_meet_test COMMAND pte_gpu_meet_test)
  set_tests_properties(pte_gpu_meet_test PROPERTIES LABELS gpu)
endif()
