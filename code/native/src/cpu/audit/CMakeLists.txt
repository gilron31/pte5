# Auditing a run at a bound nobody can cross-check: known solutions from outside
# this codebase, sampled segment re-verification against the reference path.
add_library(pte_cpu_audit STATIC audit.cc)
target_include_directories(pte_cpu_audit PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_audit
  PUBLIC pte_cpu_core pte_cpu_generate
  PRIVATE pte_cpu_meet pte_cpu_confirm pte_cpu_enumerate pte_warnings
)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_audit_test audit_test.cc)
  target_link_libraries(pte_cpu_audit_test
    PRIVATE pte_cpu_audit pte_cpu_generate pte_cpu_factorbase
            GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_audit_test COMMAND pte_cpu_audit_test)
  set_tests_properties(pte_cpu_audit_test PROPERTIES LABELS cpu)
endif()
