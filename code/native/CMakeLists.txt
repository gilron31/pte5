cmake_minimum_required(VERSION 3.24)

# LANGUAGES CXX only here; CUDA is enabled below only if a toolkit is found, so
# this project still configures on a machine with no nvcc at all.
project(pte_native LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # the python module links these static libs

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

option(PTE_ENABLE_CUDA  "Build the CUDA library, its tests and the hybrid demo" ON)
option(PTE_BUILD_TESTS  "Build the C++/CUDA unit tests" ON)
option(PTE_BUILD_PYTHON "Build the pybind11 extension module" ON)

# ---------------------------------------------------------------------------
# CUDA
# ---------------------------------------------------------------------------
# NOTE: /usr/bin/nvcc on this machine is CUDA 12.0, which hard-errors on GCC > 12.
# The presets pin /usr/local/cuda/bin/nvcc (12.9) + g++-13. See README.md.
set(PTE_HAVE_CUDA OFF)
if(PTE_ENABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES 89)  # RTX 4060 (Ada)
    endif()
    set(PTE_HAVE_CUDA ON)
    message(STATUS "CUDA enabled: ${CMAKE_CUDA_COMPILER_VERSION} "
                   "(${CMAKE_CUDA_COMPILER}), arch ${CMAKE_CUDA_ARCHITECTURES}")
  else()
    message(WARNING "No CUDA compiler found - building the CPU half only.")
  endif()
endif()

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
# Boost: CONFIG mode. Module-mode FindBoost is deprecated and this project is
# built with CMake 4. Ubuntu's libboost1.83-dev ships BoostConfig.cmake.
find_package(Boost 1.74 REQUIRED CONFIG)

# GMP has no CMake find module, but does ship pkg-config files.
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp gmpxx)

# Warning flags applied to our own targets only (never to gtest/pybind11).
add_library(pte_warnings INTERFACE)
target_compile_options(pte_warnings INTERFACE
  $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>
  $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra;--expt-relaxed-constexpr>
  $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:RelWithDebInfo,Debug>>:-lineinfo>
)

if(PTE_BUILD_TESTS)
  include(FetchContent)
  # v1.17.0: GoogleTest 1.14 declares cmake_minimum_required(3.13), which CMake 4
  # rejects outright.
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.17.0
    GIT_SHALLOW    TRUE
  )
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK   OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  enable_testing()
endif()

add_subdirectory(src/cpu)
if(PTE_HAVE_CUDA)
  add_subdirectory(src/gpu)
  add_subdirectory(apps)
endif()
if(PTE_BUILD_PYTHON)
  add_subdirectory(src/bindings)
endif()
