# The pure-C++ half: GMP + Boost, no CUDA anywhere.
add_library(pte_cpu STATIC
  bigint.cc
  matmul.cc
)
target_include_directories(pte_cpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pte_cpu
  PUBLIC  PkgConfig::GMP      # gmp + gmpxx; mpz_class appears in our headers
  PRIVATE Boost::headers pte_warnings
)

if(PTE_BUILD_TESTS)
  # One binary per library, gathering the co-located *_test.cc files. It links
  # pte_cpu and nothing else, so the CPU suite is runnable on a machine with no
  # GPU and no CUDA toolkit.
  add_executable(pte_cpu_test
    bigint_test.cc
    matmul_test.cc
  )
  target_link_libraries(pte_cpu_test
    PRIVATE pte_cpu Boost::headers GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_test COMMAND pte_cpu_test)
  set_tests_properties(pte_cpu_test PROPERTIES LABELS cpu)
endif()
