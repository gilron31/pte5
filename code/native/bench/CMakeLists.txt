# Performance benchmarks.
#
# These are *measurements*, not tests. They may be slow, they print tables
# rather than passing or failing, and nothing in src/ may depend on them. What
# they must not do is rot, so each one is registered as a ctest case at a size
# small enough to run in well under a second -- enough to catch a bench that
# stopped compiling or started crashing, not to check its numbers.
#
# Run the real thing by hand:
#   ./build/dev/bench/bench_walk_phases 1000000000 5
#
# See harness.h for why every result carries an IQR and a measured clock, and
# ceiling.h for the machine model the "%ceil" column is against.

# The harness is header-only; this target exists to carry the include path and
# the warning flags to every bench without repeating them.
add_library(pte_bench_harness INTERFACE)
target_include_directories(pte_bench_harness INTERFACE ${CMAKE_SOURCE_DIR})
target_link_libraries(pte_bench_harness INTERFACE pte_warnings)

# name;libraries;smoke-args
function(pte_add_bench name)
  cmake_parse_arguments(ARG "" "" "LIBS;SMOKE" ${ARGN})
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} PRIVATE pte_bench_harness ${ARG_LIBS})
  if(PTE_BUILD_TESTS AND ARG_SMOKE)
    add_test(NAME ${name}_smoke COMMAND ${name} ${ARG_SMOKE})
    set_tests_properties(${name}_smoke PROPERTIES LABELS bench)
  endif()
endfunction()

# B1 -- phase profile of the campaign's inner loop, by ablation.
pte_add_bench(bench_walk_phases
  LIBS pte_cpu_generate pte_cpu_factorbase pte_cpu_gaussian pte_cpu_meet
  SMOKE 200000 1)

# B2 -- point-set memory layout: AoS vs SoA vs packed. Answered (no); kept as a
# regression guard so the conclusion stays checkable.
pte_add_bench(bench_point_set_layout
  LIBS pte_cpu_generate pte_cpu_factorbase pte_cpu_gaussian
  SMOKE 200000 1)

# B3 -- the meet, decomposed. The highest-value bench: the meet is 57% of the
# inner loop at B = 1e9 and rising.
pte_add_bench(bench_meet
  LIBS pte_cpu_generate pte_cpu_factorbase pte_cpu_gaussian pte_cpu_meet
  SMOKE 200000 1)

# B4 -- the sort inside Spine::push, against a sorting-network ceiling.
pte_add_bench(bench_sort_ceiling
  LIBS pte_cpu_generate pte_cpu_factorbase pte_cpu_gaussian
  SMOKE 200000 1)

# B5 -- thread scaling, and whether the alignas(64) false-sharing fix shows up.
pte_add_bench(bench_scaling
  LIBS pte_cpu_generate pte_cpu_factorbase pte_cpu_gaussian pte_cpu_meet
       Threads::Threads
  SMOKE 200000 2)

# B6 -- factor base: Map vs Load steady-state miss rate, not startup cost.
pte_add_bench(bench_factorbase
  LIBS pte_cpu_factorbase pte_cpu_generate pte_cpu_gaussian
  SMOKE 2000000 101)

if(PTE_HAVE_CUDA)
  add_subdirectory(gpu)
endif()
