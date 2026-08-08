# The CUDA half. Each app owns a directory with its own library and test binary;
# `core` holds what they share (error checking, buffers, device discovery).
#
# Dependencies point one way only: apps depend on core, never on each other.
add_subdirectory(core)
add_subdirectory(bigint)
add_subdirectory(matmul)
add_subdirectory(meet)

# Aggregate target, so consumers can link "the GPU half" without naming apps.
add_library(pte_gpu INTERFACE)
target_link_libraries(pte_gpu INTERFACE
  pte_gpu_core pte_gpu_bigint pte_gpu_matmul pte_gpu_meet
)
