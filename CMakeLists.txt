# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.5)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

project(TERRAFERMA)
set(TERRAFERMA_VERSION_RELEASE "1")
set(TERRAFERMA_VERSION_MAJOR   "2")
set(TERRAFERMA_VERSION_MINOR   "0")
set(TERRAFERMA_VERSION_PATCH   "0")
set(TERRAFERMA_VERSION_SHORT "${TERRAFERMA_VERSION_MAJOR}.${TERRAFERMA_VERSION_MINOR}.${TERRAFERMA_VERSION_PATCH}")
set(TERRAFERMA_VERSION "${TERRAFERMA_VERSION_SHORT}")
if(NOT ${TERRAFERMA_VERSION_RELEASE})
  set(TERRAFERMA_VERSION "${TERRAFERMA_VERSION}+")
endif()

#------------------------------------------------------------------------------
# Require and use C++11

# Use C++11
set(CMAKE_CXX_STANDARD 11)

# Require C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Do not enable compler-specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

#------------------------------------------------------------------------------
# Check compiler version

# Check for GCC version - earlier versions have insuffcient C++11
# support, or bugs.
if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8.3)
    message(FATAL_ERROR "GCC version must be at least 4.8.3 (for sufficient C++11  support and to avoid some compiler bugs). You have version ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
endif()

set(BUCKETTOOLS_DIR "${PROJECT_SOURCE_DIR}/buckettools")

add_subdirectory(${BUCKETTOOLS_DIR} "buckettools")

install(FILES ${PROJECT_SOURCE_DIR}/cpp/main.cpp DESTINATION share/terraferma/cpp)
configure_file(${PROJECT_SOURCE_DIR}/cpp/CMakeLists.txt.in ${PROJECT_BINARY_DIR}/cpp/CMakeLists.txt @ONLY ESCAPE_QUOTES)
install(FILES ${PROJECT_BINARY_DIR}/cpp/CMakeLists.txt DESTINATION share/terraferma/cpp)

file(
    WRITE ${PROJECT_BINARY_DIR}/terraferma.conf
"# Helper file for setting non-default TerraFERMA and Buckettools environment variables

# Common Unix variables
export TF_CMAKE_PATH=${CMAKE_INSTALL_PREFIX}/share/terraferma/cpp
export CMAKE_MODULE_PATH=${CMAKE_INSTALL_PREFIX}/share/buckettools/cmake/modules:\$CMAKE_MODULE_PATH
export PYTHONPATH=${BUCKETTOOLS_PYTHON_INSTALL_DIR}:\$PYTHONPATH
export DIAMOND_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/buckettools/diamond:\$DIAMOND_CONFIG_PATH
export PATH=${CMAKE_INSTALL_PREFIX}/bin:\$PATH

# Special Mac variables
export DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:\$DYLD_LIBRARY_PATH
"
    )

install(FILES ${PROJECT_BINARY_DIR}/terraferma.conf DESTINATION share/terraferma)

file(
    WRITE ${PROJECT_BINARY_DIR}/terraferma.configmodule
"#%Module1.0
proc ModulesHelp { } {
        global dotversion

        puts stderr \"\tSets environment for TerraFERMA.\"
}
module-whatis \"Sets environment for TerraFERMA and Buckettools.\"

setenv       TF_CMAKE_PATH         ${CMAKE_INSTALL_PREFIX}/share/terraferma/cpp
prepend-path CMAKE_MODULE_PATH     ${CMAKE_INSTALL_PREFIX}/share/buckettools/cmake/modules
prepend-path PYTHONPATH            ${BUCKETTOOLS_PYTHON_INSTALL_DIR}
prepend-path DIAMOND_CONFIG_PATH   ${CMAKE_INSTALL_PREFIX}/share/buckettools/diamond
prepend-path PATH                  ${CMAKE_INSTALL_PREFIX}/bin

# Special Environment variables for resolving Dynamic Libraries on MacOSX
prepend-path DYLD_LIBRARY_PATH     ${CMAKE_INSTALL_PREFIX}/lib
"
    )

install(FILES ${PROJECT_BINARY_DIR}/terraferma.configmodule DESTINATION share/terraferma)

# work out some stuff used by both the tests and tutorials
find_program(HAVE_GIT git)
if(HAVE_GIT)
    execute_process (
       COMMAND git rev-parse --abbrev-ref HEAD
       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
       OUTPUT_VARIABLE GIT_BRANCH
       ERROR_VARIABLE GIT_ERR
       OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_STRIP_TRAILING_WHITESPACE
       )
    if(GIT_ERR)
       set(GIT_BRANCH "")
    endif(GIT_ERR)
else(HAVE_GIT)
    set(GIT_BRANCH "")
endif(HAVE_GIT)

set(EXCLUDE_TAGS "not:branch:${GIT_BRANCH} not:petsc:${PETSC_VERSION}")

add_subdirectory(tests)
add_subdirectory(tutorials)

# print post-install message
add_subdirectory(cmake/post-install)
