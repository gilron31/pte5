# The neutral half: types both src/cpu and src/gpu depend on, and neither owns.
#
# Header-only by construction, and the constraint that keeps it neutral is
# absolute: **no GMP, no Boost, no CUDA runtime.** A header here is included by
# translation units compiled with g++ and by translation units compiled with
# nvcc, so anything it drags in is dragged into both halves -- which is exactly
# the dependency edge REVIEW.md D5 refused for MeetHit, and the reason D5 named
# a neutral src/core/ as the right home if a shared type ever became necessary.
# S-E is that case: U256 crosses the PCIe bus as raw bytes, so the host and the
# device must agree on its layout by construction rather than by review.
add_library(pte_core INTERFACE)
target_include_directories(pte_core INTERFACE ${CMAKE_SOURCE_DIR}/src)

if(PTE_BUILD_TESTS)
  # Labelled `cpu` because it needs no device, not because it belongs to the
  # CPU half: `ctest -L cpu` means "runs on a machine with no GPU".
  #
  # Boost is linked by the *test* only, never by the library, and it is there to
  # be the independent oracle (cpp_int) in the same way cpu/bigint's test uses
  # it. The device-side agreement check for the same header lives in
  # gpu/bigint/u256_device_test.cu.
  add_executable(pte_core_test uint_test.cc window_test.cc)
  target_link_libraries(pte_core_test
    PRIVATE pte_core Boost::headers GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_core_test COMMAND pte_core_test)
  set_tests_properties(pte_core_test PROPERTIES LABELS cpu)
endif()
