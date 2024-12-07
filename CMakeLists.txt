cmake_minimum_required(VERSION 3.12)
project(kommunopp)

# Always use the Boost-provided module.
cmake_policy(SET CMP0167 NEW)

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

target_link_libraries(kommunopp PUBLIC Boost::boost)
target_link_libraries(kommunopp_test PUBLIC Boost::unit_test_framework)
