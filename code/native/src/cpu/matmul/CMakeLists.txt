add_library(pte_cpu_matmul STATIC matmul.cc)
target_include_directories(pte_cpu_matmul PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_matmul PRIVATE Boost::headers pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_matmul_test matmul_test.cc)
  target_link_libraries(pte_cpu_matmul_test
    PRIVATE pte_cpu_matmul Boost::headers GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_matmul_test COMMAND pte_cpu_matmul_test)
  set_tests_properties(pte_cpu_matmul_test PROPERTIES LABELS cpu)
endif()
