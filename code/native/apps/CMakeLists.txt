# The hybrid demo is ordinary C++ (compiled by g++, not nvcc); it reaches the
# kernels only through the plain-C++ per-app headers under src/gpu.
add_executable(hybrid_demo hybrid_demo.cc)
target_link_libraries(hybrid_demo PRIVATE pte_cpu pte_gpu pte_warnings)

# The E_4 enumeration demo: CPU generates norms and point sets, GPU does the
# meet, CPU verifies. Same seam rule -- ordinary C++, no nvcc.
add_executable(enum_demo enum_demo.cc)
target_link_libraries(enum_demo PRIVATE pte_cpu pte_gpu pte_warnings)

# The CPU-vs-GPU cross-check is the shipped correctness fixture, so it belongs
# in ctest rather than in whoever remembers to run the demo. enum_demo exits
# non-zero on any disagreement or failed confirmation.
if(PTE_BUILD_TESTS)
  add_test(NAME enum_demo_cross_check COMMAND enum_demo)
  set_tests_properties(enum_demo_cross_check PROPERTIES LABELS gpu)
endif()
