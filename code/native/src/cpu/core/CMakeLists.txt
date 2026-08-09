# Shared CPU policy: precision.h, which states the three-way exact/filter/
# certificate split every other CPU app is sized against. Header-only.
add_library(pte_cpu_core INTERFACE)
target_include_directories(pte_cpu_core INTERFACE ${CMAKE_SOURCE_DIR}/src)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_core_test precision_test.cc)
  target_link_libraries(pte_cpu_core_test
    PRIVATE pte_cpu_core GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_core_test COMMAND pte_cpu_core_test)
  set_tests_properties(pte_cpu_core_test PROPERTIES LABELS cpu)
endif()
