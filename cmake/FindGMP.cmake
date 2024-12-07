# From Z3 Sources:
# https://github.com/Z3Prover/z3/blob/master/cmake/modules/FindGMP.cmake

# License:
# Z3
# Copyright (c) Microsoft Corporation
# All rights reserved. 
# MIT License

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Tries to find an install of the GNU multiple precision library
#
# Once done this will define
#  GMP_FOUND - BOOL: System has the GMP library installed
#  GMP_INCLUDE_DIRS - LIST:The GMP include directories
#  GMP_C_LIBRARIES - LIST:The libraries needed to use GMP via it's C interface
#  GMP_CXX_LIBRARIES - LIST:The libraries needed to use GMP via it's C++ interface

include(FindPackageHandleStandardArgs)

# Try to find libraries
find_library(GMP_C_LIBRARIES
  NAMES gmp
  DOC "GMP C libraries"
)
find_library(GMP_CXX_LIBRARIES
  NAMES gmpxx
  DOC "GMP C++ libraries"
)

# Try to find headers
find_path(GMP_C_INCLUDES
  NAMES gmp.h
  DOC "GMP C header"
)

find_path(GMP_CXX_INCLUDES
  NAMES gmpxx.h
  DOC "GMP C++ header"
)

# TODO: We should check we can link some simple code against libgmp and libgmpxx

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``GMP_FOUND``
find_package_handle_standard_args(GMP
	REQUIRED_VARS GMP_C_LIBRARIES GMP_C_INCLUDES GMP_CXX_LIBRARIES GMP_CXX_INCLUDES)

if (GMP_FOUND)
  set(GMP_INCLUDE_DIRS "${GMP_C_INCLUDES}" "${GMP_CXX_INCLUDES}")
  list(REMOVE_DUPLICATES GMP_INCLUDE_DIRS)

  if (NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_C_INCLUDES}"
      IMPORTED_LOCATION "${GMP_C_LIBRARIES}")
  endif()
endif()
