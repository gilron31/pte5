# S-E, the champion meet, in the three-layer shape:
#   champion.h        plain C++ API (no CUDA types), so g++ TUs can call it
#   champion.cc       host side: sort, window planning, the session, the exact
#                     re-check of every candidate
#   champion_impl.h   device-side contract over device buffers
#   champion_impl.cu  the kernel -- the only file here that needs nvcc
add_library(pte_gpu_champion STATIC champion.cc champion_impl.cu)
target_link_libraries(pte_gpu_champion
  PUBLIC pte_gpu_core pte_core PRIVATE pte_warnings
)

if(PTE_BUILD_TESTS)
  # Goes through champion.h, links pte_gpu_champion alone -- never pte_cpu --
  # and writes its own host oracle. Being a .cc is the proof: a test that
  # cannot contain device code cannot reach past the public header.
  add_executable(pte_gpu_champion_test champion_test.cc)
  target_link_libraries(pte_gpu_champion_test
    PRIVATE pte_gpu_champion GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_gpu_champion_test COMMAND pte_gpu_champion_test)
  set_tests_properties(pte_gpu_champion_test PROPERTIES LABELS gpu)
endif()
