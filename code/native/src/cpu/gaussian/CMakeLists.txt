# Gaussian-integer arithmetic, Cornacchia, and point sets: the bottom of the
# E_4 enumeration.
add_library(pte_cpu_gaussian STATIC gaussian.cc)
target_include_directories(pte_cpu_gaussian PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_gaussian PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_gaussian_test gaussian_test.cc)
  target_link_libraries(pte_cpu_gaussian_test
    PRIVATE pte_cpu_gaussian GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_gaussian_test COMMAND pte_cpu_gaussian_test)
  set_tests_properties(pte_cpu_gaussian_test PROPERTIES LABELS cpu)
endif()
