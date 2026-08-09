# Layer 3: tower recovery, the degree-16 identity, and the divisibility ladder.
# GMP is unavoidable here -- L3 reaches 10^48 and the identity's coefficients
# ~300 bits -- and affordable, because this runs once per solution.
add_library(pte_cpu_confirm STATIC confirm.cc)
target_include_directories(pte_cpu_confirm PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_confirm
  PUBLIC pte_cpu_core pte_cpu_gaussian PkgConfig::GMP # mpz_class is in the header
  PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_confirm_test confirm_test.cc)
  target_link_libraries(pte_cpu_confirm_test
    PRIVATE pte_cpu_confirm pte_cpu_generate pte_cpu_meet pte_cpu_enumerate
            GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_confirm_test COMMAND pte_cpu_confirm_test)
  set_tests_properties(pte_cpu_confirm_test PROPERTIES LABELS cpu)
endif()
