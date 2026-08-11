# The pure-C++ half: GMP + Boost, no CUDA anywhere. Each app owns a directory
# with its own library and test binary, mirroring src/gpu.
add_subdirectory(core)
add_subdirectory(bigint)
add_subdirectory(matmul)
add_subdirectory(gaussian)
add_subdirectory(factorbase)
add_subdirectory(enumerate)
add_subdirectory(meet)
add_subdirectory(generate)
add_subdirectory(confirm)
add_subdirectory(campaign)
add_subdirectory(champion)
add_subdirectory(tower)
add_subdirectory(audit)

# Aggregate target, so consumers can link "the CPU half" without naming apps.
add_library(pte_cpu INTERFACE)
target_link_libraries(pte_cpu INTERFACE
  pte_cpu_core pte_cpu_bigint pte_cpu_matmul pte_cpu_gaussian
  pte_cpu_factorbase pte_cpu_enumerate pte_cpu_meet pte_cpu_generate
  pte_cpu_confirm pte_cpu_campaign pte_cpu_audit pte_cpu_champion
  pte_cpu_tower
)
