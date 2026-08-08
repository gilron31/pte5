add_library(pte_cpu_bigint STATIC bigint.cc)
# Include root is src/, so headers are reached as "cpu/bigint/...", keeping them
# distinct from the identically named gpu/bigint/bigint.h.
target_include_directories(pte_cpu_bigint PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_bigint
  PUBLIC  PkgConfig::GMP      # gmp + gmpxx; mpz_class appears in our header
  PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_bigint_test bigint_test.cc)
  target_link_libraries(pte_cpu_bigint_test
    PRIVATE pte_cpu_bigint Boost::headers GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_bigint_test COMMAND pte_cpu_bigint_test)
  set_tests_properties(pte_cpu_bigint_test PROPERTIES LABELS cpu)
endif()
