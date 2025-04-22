cmake_minimum_required(VERSION 3.12)
project(f4ncgb VERSION 0.1)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(ENABLE_PROFILING "Enable tracking statistics" ON)
option(ENABLE_FUZZER "Enable building the fuzzer, if available" OFF)
option(ENABLE_FAST_MERSENNE_PRIME_MODULO "Enable fast mersenne prime modulo operation" ON)
option(ENABLE_MONOMIAL_PRODUCTS_MAP "Track and cache monomial products in map (usually slower)" OFF)
option(ENABLE_COMPACT_MONOMIAL_MAP "Use a more compact map for lookups in the monomial store (usually slower)" OFF)
option(ENABLE_STORE_TRACE "Enable a (very expensive) trace of each access (operator[]) into a store." OFF)

set(MAX_BLOCKS 16 CACHE STRING "Maximum number of blocks to compile for. Must be integral. Changes compilation time (duplicates implementation for each block).")

set(F4NCGB_SRCS
  store.cpp
  signal_statistics.cpp
  parser.cpp
  monomial_trie.cpp
  primes.cpp
  linear_algebra.cpp
  zip_utils.cpp
)

set(F4NCGB_TEST_SRCS
  test.cpp
  test_simple.cpp
  test_parse.cpp
  test_fast_div.cpp
)

if(ENABLE_PROFILING)
  set(F4NCGB_SRCS ${F4NCGB_SRCS} profiling.cpp)
endif()

set(F4NCGB_TARGETS
  f4ncgb.cpp
)

if(ENABLE_STORE_TRACE)
  list(APPEND F4NCGB_SRCS store_tracer.cpp)
endif()

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

find_package(mimalloc QUIET)

# Do more optimizations in release mode. Also disable assertions.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -g -DNDEBUG")

set(CMAKE_CXX_STANDARD 20)

if(TARGET mimalloc)
  list(APPEND F4NCGB_SRCS enable-mimalloc.cpp)
  message(STATUS "Found mimalloc and use it")
else()
  message(STATUS "Not using mimalloc, but it gives you approx. 10-15% speedup. Consider installing its development package! See https://github.com/microsoft/mimalloc")
endif()

add_library(f4ncgb STATIC ${F4NCGB_SRCS})


add_executable(f4ncgb_test ${F4NCGB_TEST_SRCS})
add_executable(f4ncgb_bin ${F4NCGB_TARGETS})
target_link_libraries(f4ncgb_test PUBLIC f4ncgb)
target_link_libraries(f4ncgb_bin PUBLIC f4ncgb)
set_target_properties(f4ncgb_bin PROPERTIES OUTPUT_NAME "f4ncgb")

set_property(TARGET f4ncgb PROPERTY CXX_STANDARD 20)
set_property(TARGET f4ncgb_test PROPERTY CXX_STANDARD 20)
set_property(TARGET f4ncgb_bin PROPERTY CXX_STANDARD 20)

if(TARGET mimalloc)
  target_link_libraries(f4ncgb PUBLIC mimalloc)
endif()

target_compile_definitions(f4ncgb PUBLIC F4NCGB_MAX_BLOCKS=${MAX_BLOCKS})
message(STATUS "Compiling with " ${MAX_BLOCKS} " maximum blocks.")

if(ENABLE_FAST_MERSENNE_PRIME_MODULO)
  message(STATUS "Fast % 2^31-1 (mersenne prime) modulo operation is enabled.")
  target_compile_definitions(f4ncgb PUBLIC F4NCGB_FAST_MERSENNE_PRIME_MODULO)
else()
  message(STATUS "Fast % 2^31-1 (mersenne prime) modulo operation is disabled.")
endif()

if(ENABLE_MONOMIAL_PRODUCTS_MAP)
  message(STATUS "Monomial products map is enabled.")
  target_compile_definitions(f4ncgb PUBLIC F4NCGB_USE_MONOMIAL_PRODUCTS_MAP)
else()
  message(STATUS "Monomial products map is disabled.")
endif()

if(ENABLE_COMPACT_MONOMIAL_MAP)
  message(STATUS "Compact monomial map is enabled")
  target_compile_definitions(f4ncgb PUBLIC F4NCGB_USE_COMPACT_MONOMIAL_MAP)
else()
  message(STATUS "Compact monomial map is disabled")
endif()

if(ENABLE_STORE_TRACE)
  message(STATUS "Store tracing is enabled")
  target_compile_definitions(f4ncgb PUBLIC F4NCGB_ENABLE_STORE_TRACE)
else()
  message(STATUS "Store tracing is disabled")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Sanitize")
  target_compile_options(f4ncgb PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(f4ncgb PRIVATE -fsanitize=address -fsanitize=undefined)
  target_compile_options(f4ncgb_bin PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(f4ncgb_bin PRIVATE -fsanitize=address -fsanitize=undefined)
  target_compile_options(f4ncgb_test PRIVATE -fsanitize=address -fsanitize=undefined -g)
  target_link_options(f4ncgb_test PRIVATE -fsanitize=address -fsanitize=undefined)

  target_link_libraries(f4ncgb_test PRIVATE -fsanitize=address -fsanitize=undefined)
  target_link_libraries(f4ncgb_bin PRIVATE -fsanitize=address -fsanitize=undefined)
endif()

find_package(Boost COMPONENTS unit_test_framework program_options REQUIRED)
find_package(GMP REQUIRED)
find_package(Flint REQUIRED)

target_compile_options(f4ncgb PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_compile_options(f4ncgb_test PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_compile_options(f4ncgb_bin PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wall -Wextra -Wsign-conversion>
     $<$<CXX_COMPILER_ID:MSVC>:
          /W4>)

target_link_libraries(f4ncgb PUBLIC Boost::boost GMP::GMP)
target_link_libraries(f4ncgb_test PUBLIC Boost::unit_test_framework)
target_link_libraries(f4ncgb PUBLIC Flint::Flint Boost::boost GMP::GMP)
target_link_libraries(f4ncgb_bin PUBLIC Boost::boost Boost::program_options)

target_include_directories(f4ncgb PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(f4ncgb_bin PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_PROFILING)
  target_compile_definitions(f4ncgb PUBLIC F4NCGB_ENABLE_PROFILING)

  # This target makes no sense otherwise.
  add_executable(mod_prime_tester mod_prime_tester.cpp)
  target_link_libraries(mod_prime_tester PUBLIC Boost::boost Boost::program_options)
  target_include_directories(mod_prime_tester PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  target_compile_definitions(mod_prime_tester PUBLIC F4NCGB_ENABLE_PROFILING)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported AND CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "IPO / LTO enabled (Release build and compiler supports it)")
  set_property(TARGET f4ncgb_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fuzztest/CMakeLists.txt" AND ENABLE_FUZZER)
  # Create the fuzzing target only if the fuzztest subdirectory
  # exists.
  message(STATUS "Fuzzer enabled.")
  set(FUZZTEST_FUZZING_MODE ON)
  add_subdirectory(fuzztest)
  fuzztest_setup_fuzzing_flags()

  add_executable(
    f4ncgb_fuzz
    fuzzer.cpp
  )

  link_fuzztest(f4ncgb_fuzz)
  target_link_libraries(f4ncgb_fuzz PUBLIC f4ncgb)
  target_compile_options(f4ncgb_fuzz PRIVATE -g)

  if(CMAKE_BUILD_TYPE MATCHES "Sanitize")
    target_link_options(f4ncgb_fuzz PRIVATE -fsanitize=address -fsanitize=undefined)
    target_link_libraries(f4ncgb_fuzz PRIVATE -fsanitize=address -fsanitize=undefined)
  endif()
else()
  message(STATUS "Fuzzer disabled.")
endif()
