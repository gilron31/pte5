# The split-prime factor base: segmented sieve, decomposition, artifact I/O.
# Segmented so resident memory is constant (0.2 MB) at any bound -- see the
# header for the measured table.
add_library(pte_cpu_factorbase STATIC factorbase.cc)
target_include_directories(pte_cpu_factorbase PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_factorbase
  PUBLIC pte_cpu_core pte_cpu_gaussian PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_factorbase_test factorbase_test.cc)
  # Links pte_cpu_enumerate only for split_primes_up_to, the monolithic sieve
  # this one replaces: agreement with it is the fixture that makes the swap safe.
  target_link_libraries(pte_cpu_factorbase_test
    PRIVATE pte_cpu_factorbase pte_cpu_enumerate GTest::gtest_main pte_warnings
  )
  # So the test can find the shipped artifact without hardcoding a path.
  target_compile_definitions(pte_cpu_factorbase_test
    PRIVATE PTE_DATA_DIR="${CMAKE_SOURCE_DIR}/../../data"
  )
  add_test(NAME pte_cpu_factorbase_test COMMAND pte_cpu_factorbase_test)
  set_tests_properties(pte_cpu_factorbase_test PROPERTIES LABELS cpu)
endif()
