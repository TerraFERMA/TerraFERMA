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

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)

# Default build type (can be overridden by user)
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
  "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo.")

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Add Mac specific loader flags
if (APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,dynamic_lookup" CACHE STRING "should only be set if macosx" FORCE)
endif()

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

project(terraferma)

set(EXECUTABLE "terraferma" CACHE STRING
     "The name of the executable")

set(LOGLEVEL "INFO" CACHE STRING
     "The log level used in the run target, options are: DBG, TRACE, PROGRESS, INFO, WARNING, ERROR, CRITICAL.")

set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})

find_package(BUCKETTOOLS REQUIRED)

if (NOT BUCKETTOOLS_FOUND)
  message(FATAL_ERROR "Could not find BucketTools.")
endif()

include_directories(${BUCKETTOOLS_INCLUDE_DIRS})

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
                           "import sys; sys.stdout.write(sys.version[:5])"
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

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

set(OPTIONSFILE CACHE FILEPATH "Choose the input file.")

if(NOT IS_ABSOLUTE ${OPTIONSFILE})
  set(OPTIONSFILE ${PROJECT_BINARY_DIR}/${OPTIONSFILE})
endif(NOT IS_ABSOLUTE ${OPTIONSFILE})

if (IS_DIRECTORY ${OPTIONSFILE})
  message(FATAL_ERROR "OPTIONSFILE not set or set to a directory.")
endif(IS_DIRECTORY ${OPTIONSFILE})

if (NOT EXISTS ${OPTIONSFILE})
  message(FATAL_ERROR "Could not find OPTIONSFILE ${OPTIONSFILE}")
endif(NOT EXISTS ${OPTIONSFILE})

add_subdirectory(${BUCKETTOOLS_SHARED_DIRS}/buckettools/ufc "buckettools_ufc")

add_executable(${EXECUTABLE} ${PROJECT_SOURCE_DIR}/main.cpp)

set (EXTRA_LIBS ${EXTRA_LIBS} ${BUCKETTOOLS_LIBRARIES} buckettools_ufc)

target_link_libraries(${EXECUTABLE} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${EXTRA_LIBS} ${SPUD_LIBRARIES} ${PYTHON_LIBRARIES})

add_custom_target(
  run
  COMMAND ${EXECUTABLE} -v${LOGLEVEL} -l ${OPTIONSFILE} 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${EXECUTABLE}
  )

