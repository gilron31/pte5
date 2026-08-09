# The production meet: compile-time partition tables for D <= 8, an O(D^2)
# hashed path for the large-D tail that carries 74% of the work. Header-only --
# every entry point is a template on the key width and the sink, so that one
# source can serve the host and (later) the device.
add_library(pte_cpu_meet INTERFACE)
target_include_directories(pte_cpu_meet INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_meet INTERFACE pte_cpu_core)

if(PTE_BUILD_TESTS)
  # Links pte_cpu_enumerate for meet_oracle and the Enumerator. That is the
  # dependency the cross-check is made of, and the only place the two meets are
  # allowed to touch: the oracle must stay independently written.
  add_executable(pte_cpu_meet_test meet_test.cc)
  target_link_libraries(pte_cpu_meet_test
    PRIVATE pte_cpu_meet pte_cpu_enumerate GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_meet_test COMMAND pte_cpu_meet_test)
  set_tests_properties(pte_cpu_meet_test PROPERTIES LABELS cpu)
endif()
