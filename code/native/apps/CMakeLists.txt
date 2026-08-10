# The exhaustive E_4 campaign. Pure CPU -- no CUDA anywhere -- so it builds and
# runs on a machine with no device, which is the point: the campaign is CPU
# work (drafts/implementation_design.md 3.1) and the GPU's value in this project
# is the S-D and S-E searches.
add_executable(enum_campaign enum_campaign.cc)
target_link_libraries(enum_campaign PRIVATE pte_cpu_campaign pte_cpu_audit pte_warnings)

if(PTE_HAVE_CUDA)
  # The hybrid demo is ordinary C++ (compiled by g++, not nvcc); it reaches the
  # kernels only through the plain-C++ per-app headers under src/gpu.
  add_executable(hybrid_demo hybrid_demo.cc)
  target_link_libraries(hybrid_demo PRIVATE pte_cpu pte_gpu pte_warnings)

  # The E_4 enumeration demo: CPU generates norms and point sets, GPU does the
  # meet, CPU verifies and certifies.
  add_executable(enum_demo enum_demo.cc)
  target_link_libraries(enum_demo PRIVATE pte_cpu pte_gpu pte_warnings)

  # S-E, the champion meet: one enormous point set, the device hash-insert
  # kernel, and the CPU reference checking it. Takes the prime count as an
  # argument -- 20 is the sized champion (D = 2^19), the default 14 runs in
  # seconds.
  add_executable(champion_demo champion_demo.cc)
  target_link_libraries(champion_demo
    PRIVATE pte_cpu_champion pte_gpu_champion pte_gpu_core pte_warnings
  )

  # The CPU-vs-GPU cross-checks are the shipped correctness fixtures, so they
  # belong in ctest rather than in whoever remembers to run the demos. Both
  # exit non-zero on any disagreement, failed confirmation, or failed
  # certificate.
  if(PTE_BUILD_TESTS)
    add_test(NAME enum_demo_cross_check COMMAND enum_demo)
    set_tests_properties(enum_demo_cross_check PROPERTIES LABELS gpu)
    # 12 primes = D = 2048 = 2.1e6 pairs, ~1 s including the CUDA context.
    add_test(NAME champion_demo_cross_check COMMAND champion_demo 12)
    set_tests_properties(champion_demo_cross_check PROPERTIES LABELS gpu)
  endif()
endif()

# A short end-to-end run of the shipped binary, so the CLI, the checkpoint and
# the resume path are exercised as a program rather than only as a library.
if(PTE_BUILD_TESTS)
  add_test(NAME enum_campaign_smoke
           COMMAND enum_campaign --bound=1000000 --threads=2 --batch=64 --quiet
                   --out=${CMAKE_CURRENT_BINARY_DIR}/smoke_catalog.csv)
  set_tests_properties(enum_campaign_smoke PROPERTIES LABELS cpu)
endif()
