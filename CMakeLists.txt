cmake_minimum_required(VERSION 3.12)
project(kommunopp)

# Always use the Boost-provided module.
cmake_policy(SET CMP0167 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(KOMMUNOPP_SRCS
  kommunopp.cpp
)

set(KOMMUNOPP_TEST_SRCS
  test.cpp
  test_simple.cpp
)

add_library(kommunopp STATIC ${KOMMUNOPP_SRCS})

add_executable(kommunopp_test ${KOMMUNOPP_TEST_SRCS})
target_link_libraries(kommunopp_test PUBLIC kommunopp)

set_property(TARGET kommunopp PROPERTY CXX_STANDARD 20)
set_property(TARGET kommunopp_test PROPERTY CXX_STANDARD 20)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
find_package(GMP REQUIRED)

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

target_link_libraries(kommunopp PUBLIC Boost::boost GMP::GMP)
target_link_libraries(kommunopp_test PUBLIC Boost::unit_test_framework)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "IPO / LTO enabled (Release build and compiler supports it)")
    set_property(TARGET kommunopp_test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()
