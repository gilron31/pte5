# The hybrid demo is ordinary C++ (compiled by g++, not nvcc); it reaches the
# kernels only through src/gpu/gpu.h.
add_executable(hybrid_demo hybrid_demo.cc)
target_link_libraries(hybrid_demo PRIVATE pte_cpu pte_gpu pte_warnings)
