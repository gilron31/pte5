# The hybrid demo is ordinary C++ (compiled by g++, not nvcc); it reaches the
# kernels only through src/gpu/gpu.h.
add_executable(hybrid_demo hybrid_demo.cc)
target_link_libraries(hybrid_demo PRIVATE pte_cpu pte_gpu pte_warnings)

# The E_4 enumeration demo: CPU generates norms and point sets, GPU does the
# meet, CPU verifies. Same seam rule -- ordinary C++, no nvcc.
add_executable(enum_demo enum_demo.cc)
target_link_libraries(enum_demo PRIVATE pte_cpu pte_gpu pte_warnings)
