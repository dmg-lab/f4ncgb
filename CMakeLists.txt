cmake_minimum_required(VERSION 3.12)
project(kommunopp)

# Always use the Boost-provided module.
cmake_policy(SET CMP0167 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(KOMMUNOPP_SRCS
  kommunopp.cpp
  signal_statistics.cpp
  parser.cpp
)

set(KOMMUNOPP_TEST_SRCS
  test.cpp
  test_simple.cpp
  test_parse.cpp
)

set(KOMMUNOPP_TARGETS
  freegb.cpp
)

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


set(CMAKE_C_FLAGS_SANITIZE
    "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" CACHE STRING
    "Flags used by the C compiler for Sanitize build type or configuration." FORCE)

set(CMAKE_CXX_FLAGS_SANITIZE
    "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined" CACHE STRING
    "Flags used by the C++ compiler for Sanitize build type or configuration." FORCE)

set(CMAKE_EXE_LINKER_FLAGS_SANITIZE
    "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined" CACHE STRING
    "Linker flags to be used to create executables for Sanitize build type." FORCE)

set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE
    "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined" CACHE STRING
    "Linker lags to be used to create shared libraries for Sanitize build type." FORCE)

add_library(kommunopp STATIC ${KOMMUNOPP_SRCS})

add_executable(kommunopp_test ${KOMMUNOPP_TEST_SRCS})
add_executable(freegb ${KOMMUNOPP_TARGETS})
target_link_libraries(kommunopp_test PUBLIC kommunopp)
target_link_libraries(freegb PUBLIC kommunopp)


set_property(TARGET kommunopp PROPERTY CXX_STANDARD 20)
set_property(TARGET kommunopp_test PROPERTY CXX_STANDARD 20)
set_property(TARGET freegb PROPERTY CXX_STANDARD 20)


find_package(Boost COMPONENTS unit_test_framework REQUIRED)
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
target_link_libraries(freegb PUBLIC Boost::boost GMP::GMP)
target_link_libraries(freegb PUBLIC Flint::Flint)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported AND CMAKE_BUILD_TYPE MATCHES "Release")
    message(STATUS "IPO / LTO enabled (Release build and compiler supports it)")
    set_property(TARGET kommunopp_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()
