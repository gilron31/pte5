# S-E, the champion meet -- the CPU reference.
#
# Depends on core (the fixed-width type and the value-window partition) and on
# nothing else: no GMP, no factor base, no Enumerator. That is deliberate --
# this is one function over one array of nu values, and keeping it that way is
# what lets the *test* link the enumerator and the production meet without the
# library itself acquiring them.
add_library(pte_cpu_champion STATIC champion.cc champion_points.cc)
target_include_directories(pte_cpu_champion PUBLIC ${CMAKE_SOURCE_DIR}/src)
# champion.cc needs core alone; champion_points.cc needs GMP and the factor-base
# decomposition, because a champion's coordinates are 61 bits and the int32
# point-set path would overflow silently.
target_link_libraries(pte_cpu_champion
  PUBLIC pte_core pte_cpu_bigint pte_cpu_gaussian PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  # Links the enumerator and the production meet because the fixture *is* the
  # cross-check against them: three independent implementations of one
  # question, agreeing on real point sets.
  add_executable(pte_cpu_champion_test champion_test.cc champion_points_test.cc)
  target_link_libraries(pte_cpu_champion_test
    PRIVATE pte_cpu_champion pte_cpu_enumerate pte_cpu_meet pte_cpu_gaussian
            GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_champion_test COMMAND pte_cpu_champion_test)
  set_tests_properties(pte_cpu_champion_test PROPERTIES LABELS cpu)
endif()
