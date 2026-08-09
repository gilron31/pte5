# The pure-C++ half: GMP + Boost, no CUDA anywhere. Each app owns a directory
# with its own library and test binary, mirroring src/gpu.
add_subdirectory(core)
add_subdirectory(bigint)
add_subdirectory(matmul)
add_subdirectory(gaussian)
add_subdirectory(enumerate)

# Aggregate target, so consumers can link "the CPU half" without naming apps.
add_library(pte_cpu INTERFACE)
target_link_libraries(pte_cpu INTERFACE
  pte_cpu_core pte_cpu_bigint pte_cpu_matmul pte_cpu_gaussian pte_cpu_enumerate
)
