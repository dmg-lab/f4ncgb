cmake_minimum_required(VERSION 3.12)
project(kommunopp VERSION 0.1)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(ENABLE_PROFILING "Enable tracking statistics" ON)
option(ENABLE_FUZZER "Enable building the fuzzer, if available" ON)
option(ENABLE_MONOMIAL_PRODUCTS_MAP "Track and cache monomial products in map" OFF)

set(KOMMUNOPP_SRCS
  kommunopp.cpp
  signal_statistics.cpp
  parser.cpp
  monomial_trie.cpp
  primes.cpp
  linear_algebra.cpp
)

set(KOMMUNOPP_TEST_SRCS
  test.cpp
  test_simple.cpp
  test_parse.cpp
)

if(ENABLE_PROFILING)
  set(KOMMUNOPP_SRCS ${KOMMUNOPP_SRCS} profiling.cpp)
endif()

set(KOMMUNOPP_TARGETS
  freegb.cpp
)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config.hpp")

if(isMultiConfig)
    if(NOT "Sanitize" IN_LIST CMAKE_CONFIGURATION_TYPES)
        list(APPEND CMAKE_CONFIGURATION_TYPES Sanitize)
    endif()
else()
    set(allowedBuildTypes Sanitize Debug Release RelWithDebInfo MinSizeRel)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${allowedBuildTypes}")

    if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
        message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}")
    endif()
endif()

# Do more optimizations in release mode.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

set(CMAKE_CXX_STANDARD 20)

add_library(kommunopp STATIC ${KOMMUNOPP_SRCS})

add_executable(kommunopp_test ${KOMMUNOPP_TEST_SRCS})
add_executable(freegb ${KOMMUNOPP_TARGETS})
target_link_libraries(kommunopp_test PUBLIC kommunopp)
target_link_libraries(freegb PUBLIC kommunopp)

set_property(TARGET kommunopp PROPERTY CXX_STANDARD 20)
set_property(TARGET kommunopp_test PROPERTY CXX_STANDARD 20)
set_property(TARGET freegb PROPERTY CXX_STANDARD 20)

if(ENABLE_MONOMIAL_PRODUCTS_MAP)
  message(STATUS "Monomial products map is enabled.")
  target_compile_definitions(kommunopp PUBLIC KOMMUNOPP_USE_MONOMIAL_PRODUCTS_MAP)
else()
  message(STATUS "Monomial products map is disabled.")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Sanitize")
  target_compile_options(kommunopp PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(kommunopp PRIVATE -fsanitize=address -fsanitize=undefined)
  target_compile_options(freegb PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(freegb PRIVATE -fsanitize=address -fsanitize=undefined)
  target_compile_options(kommunopp_test PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(kommunopp_test PRIVATE -fsanitize=address -fsanitize=undefined)

  target_link_libraries(kommunopp_test PRIVATE -fsanitize=address -fsanitize=undefined)
  target_link_libraries(freegb PRIVATE -fsanitize=address -fsanitize=undefined)
endif()

find_package(Boost COMPONENTS unit_test_framework program_options REQUIRED)
find_package(GMP REQUIRED)
find_package(Flint REQUIRED)

target_compile_options(kommunopp PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_compile_options(kommunopp_test PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_compile_options(freegb PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_link_libraries(kommunopp PUBLIC Boost::boost GMP::GMP)
target_link_libraries(kommunopp_test PUBLIC Boost::unit_test_framework)
target_link_libraries(kommunopp PUBLIC Flint::Flint Boost::boost GMP::GMP)
target_link_libraries(freegb PUBLIC Boost::boost Boost::program_options)

target_include_directories(freegb PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_PROFILING)
  target_compile_definitions(kommunopp PUBLIC KOMMUNOPP_ENABLE_PROFILING)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported AND CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "IPO / LTO enabled (Release build and compiler supports it)")
  set_property(TARGET kommunopp_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fuzztest/CMakeLists.txt")
  # Create the fuzzing target only if the fuzztest subdirectory
  # exists.
  set(FUZZTEST_FUZZING_MODE ON)
  add_subdirectory(fuzztest)
  fuzztest_setup_fuzzing_flags()

  add_executable(
    kommunopp_fuzz
    fuzzer.cpp
  )

  link_fuzztest(kommunopp_fuzz)
  target_link_libraries(kommunopp_fuzz PUBLIC kommunopp)
  target_compile_options(kommunopp_fuzz PRIVATE -g)

  if(CMAKE_BUILD_TYPE MATCHES "Sanitize")
    target_link_options(kommunopp_fuzz PRIVATE -fsanitize=address -fsanitize=undefined)
    target_link_libraries(kommunopp_fuzz PRIVATE -fsanitize=address -fsanitize=undefined)
  endif()
endif()
