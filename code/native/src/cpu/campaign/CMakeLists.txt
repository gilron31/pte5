# The exhaustive campaign: segments, threads, atomic checkpoints, streaming
# catalog. Links the whole CPU pipeline because it is the thing that runs it.
add_library(pte_cpu_campaign STATIC campaign.cc)
target_include_directories(pte_cpu_campaign PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pte_cpu_campaign
  PUBLIC pte_cpu_core pte_cpu_generate pte_cpu_factorbase
  PRIVATE pte_cpu_meet pte_cpu_confirm pte_cpu_enumerate Threads::Threads pte_warnings
)

if(PTE_BUILD_TESTS)
  add_executable(pte_cpu_campaign_test campaign_test.cc)
  target_link_libraries(pte_cpu_campaign_test
    PRIVATE pte_cpu_campaign pte_cpu_generate pte_cpu_meet pte_cpu_confirm
            pte_cpu_enumerate GTest::gtest_main pte_warnings
  )
  add_test(NAME pte_cpu_campaign_test COMMAND pte_cpu_campaign_test)
  set_tests_properties(pte_cpu_campaign_test PROPERTIES LABELS cpu)
endif()
