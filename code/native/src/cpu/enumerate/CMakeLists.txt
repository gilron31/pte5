# Norm generation, the reference meet, and the primitivity rule.
add_library(pte_cpu_enumerate STATIC enumerate.cc)
target_include_directories(pte_cpu_enumerate PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_enumerate PUBLIC pte_cpu_gaussian PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_enumerate_test enumerate_test.cc)
  target_link_libraries(pte_cpu_enumerate_test
    PRIVATE pte_cpu_enumerate GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_enumerate_test COMMAND pte_cpu_enumerate_test)
  set_tests_properties(pte_cpu_enumerate_test PROPERTIES LABELS cpu)
endif()
