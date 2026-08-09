# Norm generation with point sets carried down the DFS spine, replacing the
# per-norm rebuild. Fixed capacity: one allocation per Spine, none per norm.
add_library(pte_cpu_generate STATIC generate.cc)
target_include_directories(pte_cpu_generate PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_generate
  PUBLIC pte_cpu_core pte_cpu_gaussian pte_cpu_factorbase PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  # Links pte_cpu_enumerate for the reference enumerator, which rebuilds every
  # point set from scratch. That cross-check is the point of the test; the two
  # must stay independently written.
  add_executable(pte_cpu_generate_test generate_test.cc)
  target_link_libraries(pte_cpu_generate_test
    PRIVATE pte_cpu_generate pte_cpu_enumerate GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_generate_test COMMAND pte_cpu_generate_test)
  set_tests_properties(pte_cpu_generate_test PROPERTIES LABELS cpu)
endif()
