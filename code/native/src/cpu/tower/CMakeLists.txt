# G1, the local model: the depth-d tower recursion over Z/p^k.
#
# Links nothing. Not "links little" -- nothing at all beyond the include path:
# every value here is a residue mod p^k <= 2^20 and every state is a bitmask, so
# there is no GMP, no Boost and no core type in sight. That is worth preserving,
# because the same DP is the candidate S-D device kernel
# (drafts/gpu_implementation_plan.md 2.3) and a dependency acquired here would
# have to be paid for again on the device.
add_library(pte_cpu_tower STATIC tower.cc tower_oracle.cc)
target_include_directories(pte_cpu_tower PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_tower PRIVATE pte_warnings)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_tower_test tower_test.cc)
  target_link_libraries(pte_cpu_tower_test
    PRIVATE pte_cpu_tower GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_tower_test COMMAND pte_cpu_tower_test)
  set_tests_properties(pte_cpu_tower_test PROPERTIES LABELS cpu)
endif()
