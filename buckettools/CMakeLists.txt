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

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8) 

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# The project name for buckettools (sets up environment variables for binary and source directories)
project(buckettools)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)

# Need to get VTK config because VTK uses advanced VTK features which
# mean it's not enough to just link to the DOLFIN target. See
# http://www.vtk.org/pipermail/vtk-developers/2013-October/014402.html
find_package(VTK HINTS ${VTK_DIR} $ENV{VTK_DIR} NO_MODULE QUIET)

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Compiler flags
set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

# Add Mac specific loader flags
if (APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,dynamic_lookup" CACHE STRING "should only be set if macosx" FORCE)
endif()

find_package(PythonInterp 2 REQUIRED)

# Set variables to help find Python library that is compatible with
# interpreter
if (PYTHONINTERP_FOUND)
  # Get Python include path from Python interpretter
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                         "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)

  # Get Python library path from interpreter
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                        "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
                  OUTPUT_VARIABLE _PYTHON_LIB_PATH
                  RESULT_VARIABLE _PYTHON_LIB_RESULT)

  # Get Python version from interpreter (find_package sets this but only for cmake >= 2.8.5, which we don't require)
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                           "import sys; sys.stdout.write(sys.version.split()[0])"
                  OUTPUT_VARIABLE _PYTHON_VERSION_STR
                  RESULT_VARIABLE _PYTHON_VERSION_RESULT)

  # Set include path, if returned by interpreter
  if ("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
    set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH})
  else()
    message(FATAL_ERROR "Failed to find python include path.")
  endif()

  # Add a search path for Python library based on output from
  # iterpreter
  if ("${_PYTHON_LIB_RESULT}" STREQUAL "0")
    set(CMAKE_LIBRARY_PATH ${_PYTHON_LIB_PATH})
  else()
    message(FATAL_ERROR "Failed to find python library path.")
  endif()

  # Add version number based on output from iterpreter
  if ("${_PYTHON_VERSION_RESULT}" STREQUAL "0")
    set(PYTHON_VERSION_STR ${_PYTHON_VERSION_STR})
  else()
    message(FATAL_ERROR "Failed to find python version.")
  endif()
else()
  message(FATAL_ERROR "Could not find python interpreter.")
endif()

# Find Pythons livs
find_package(PythonLibs ${PYTHON_VERSION_STR} EXACT REQUIRED)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

add_subdirectory(cpp)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

if (PYTHONINTERP_FOUND)

  set(PYTHON_LIB_DISCOVER_STR "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))")

  # Probe Python interpreter
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${PYTHON_LIB_DISCOVER_STR}"
    RESULT_VARIABLE _BUCKETTOOLS_PYTHON_INSTALL_RESULT
    OUTPUT_VARIABLE _BUCKETTOOLS_PYTHON_INSTALL_DIR
    )
  
  if ("${_BUCKETTOOLS_PYTHON_INSTALL_RESULT}" STREQUAL "0")
    set(BUCKETTOOLS_PYTHON_INSTALL_DIR ${_BUCKETTOOLS_PYTHON_INSTALL_DIR} PARENT_SCOPE)
    set(BUCKETTOOLS_PYTHON_INSTALL_DIR ${_BUCKETTOOLS_PYTHON_INSTALL_DIR})
  else()
    message(FATAL_ERROR "Failed to find python installation path.")
  endif()
endif()

if (BUCKETTOOLS_PYTHON_INSTALL_DIR)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/buckettools DESTINATION ${BUCKETTOOLS_PYTHON_INSTALL_DIR})
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/schemas DESTINATION share/buckettools)

# Write a schemata file for tfmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/tfml
    "TerraFERMA Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/terraferma.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/tfml DESTINATION share/buckettools/diamond/schemata)

# Write a schemata file for tfmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/shml
    "TerraFERMA Simulation Harness Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/simulations.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/shml DESTINATION share/buckettools/diamond/schemata)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake DESTINATION share/buckettools)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/ufc   DESTINATION share/buckettools)

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/scripts/*)
file(GLOB EXCLUDESCRIPTS ${PROJECT_SOURCE_DIR}/scripts/.*)
if(EXCLUDESCRIPTS)
  list(REMOVE_ITEM SCRIPTS ${EXCLUDESCRIPTS})
endif()
install(FILES ${SCRIPTS} DESTINATION bin 
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

configure_file(${PROJECT_SOURCE_DIR}/cmake/templates/BUCKETTOOLSConfig.cmake.in
  ${PROJECT_BINARY_DIR}/BUCKETTOOLSConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/BUCKETTOOLSConfig.cmake DESTINATION share/buckettools/cmake)

