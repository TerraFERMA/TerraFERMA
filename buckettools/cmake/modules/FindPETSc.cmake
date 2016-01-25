# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDE_DIRS - include directories for PETSc
#  PETSC_LIBRARIES    - libraries for PETSc
#  PETSC_DIR          - directory where PETSc is built
#  PETSC_ARCH         - architecture for which PETSc is built
#  PETSC_VERSION      - version for PETSc
#  PETSC_VERSION_MAJOR - First number in PETSC_VERSION
#  PETSC_VERSION_MINOR - Second number in PETSC_VERSION
#  PETSC_VERSION_SUBMINOR - Third number in PETSC_VERSION
#
# This config script is (very loosley) based on a PETSc CMake script by Jed Brown.

#=============================================================================
# Copyright (C) 2010-2012 Garth N. Wells, Anders Logg and Johannes Ring
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# NOTE: The PETSc Makefile returns a bunch of libraries with '-L' and '-l',
# wheres we would prefer complete paths. For a discussion, see
# http://www.cmake.org/Wiki/CMake:Improving_Find*_Modules#Current_workarounds

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

# This file is a trimmed down version of FindPETSc taken from DOLFIN (see upper copyright notice).

message(STATUS "Checking for package 'PETSc'")

# Set debian_arches (PETSC_ARCH for Debian-style installations)
foreach (debian_arches linux kfreebsd)
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(DEBIAN_FLAVORS ${debian_arches}-gnu-c-debug ${debian_arches}-gnu-c-opt ${DEBIAN_FLAVORS})
  else()
    set(DEBIAN_FLAVORS ${debian_arches}-gnu-c-opt ${debian_arches}-gnu-c-debug ${DEBIAN_FLAVORS})
  endif()
endforeach()

# List of possible locations for PETSC_DIR
set(petsc_dir_locations "")
list(APPEND petsc_dir_locations "/usr/lib/petscdir/3.4.2")    # Debian location
list(APPEND petsc_dir_locations "/usr/lib/petscdir/3.2")    # Debian location
list(APPEND petsc_dir_locations "/usr/lib/petscdir/3.1")    # Debian location
list(APPEND petsc_dir_locations "/usr/lib/petscdir/3.0.0")  # Debian location
list(APPEND petsc_dir_locations "/opt/local/lib/petsc")     # Macports location
list(APPEND petsc_dir_locations "/usr/local/lib/petsc")     # User location
list(APPEND petsc_dir_locations "$ENV{HOME}/petsc")         # User location

# Add other possible locations for PETSC_DIR
set(_SYSTEM_LIB_PATHS "${CMAKE_SYSTEM_LIBRARY_PATH};${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
string(REGEX REPLACE ":" ";" libdirs ${_SYSTEM_LIB_PATHS})
foreach (libdir ${libdirs})
  get_filename_component(petsc_dir_location "${libdir}/" PATH)
  list(APPEND petsc_dir_locations ${petsc_dir_location})
endforeach()

# Try to figure out PETSC_DIR by finding petsc.h
find_path(PETSC_DIR include/petsc.h
  HINTS ${PETSC_DIR} $ENV{PETSC_DIR}
  PATHS ${petsc_dir_locations}
  DOC "PETSc directory")

# Report result of search for PETSC_DIR
if (DEFINED PETSC_DIR)
  message(STATUS "PETSC_DIR is ${PETSC_DIR}")
else()
  message(STATUS "PETSC_DIR is empty")
endif()

# Try to figure out PETSC_ARCH if not set
if (PETSC_DIR AND NOT PETSC_ARCH)
  set(_petsc_arches
    $ENV{PETSC_ARCH}   # If set, use environment variable first
    ${DEBIAN_FLAVORS}  # Debian defaults
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
  set(petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_petsc_arches})
    if (NOT PETSC_ARCH)
      find_path(petscconf petscconf.h
      HINTS ${PETSC_DIR}
      PATH_SUFFIXES ${arch}/include bmake/${arch}
      NO_DEFAULT_PATH)
      if (petscconf)
        set(PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture")
      endif()
    endif()
  endforeach()
  set(petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif()

# Report result of search for PETSC_ARCH
if (DEFINED PETSC_ARCH)
  message(STATUS "PETSC_ARCH is ${PETSC_ARCH}")
else()
  message(STATUS "PETSC_ARCH is empty")
endif()

# Look for petscconf.h
if (EXISTS ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h)
  message(STATUS "Found petscconf.h")
  set(FOUND_PETSC_CONF 1)
else()
  message(STATUS "Unable to find petscconf.h")
endif()

# Get variables from PETSc configuration
if (FOUND_PETSC_CONF)

  # Find PETSc config file
  find_file(PETSC_VARIABLES_FILE NAMES variables PATHS ${PETSC_DIR}/lib/petsc/conf
    ${PETSC_DIR}/lib/petsc-conf ${PETSC_DIR}/conf)

  # Create a temporary Makefile to probe the PETSc configuration
  set(petsc_config_makefile ${PROJECT_BINARY_DIR}/Makefile.petsc)
  file(WRITE ${petsc_config_makefile}
"# This file was autogenerated by FindPETSc.cmake
PETSC_DIR  = ${PETSC_DIR}
PETSC_ARCH = ${PETSC_ARCH}
-include ${PETSC_VARIABLES_FILE}
show :
	-@echo -n \${\${VARIABLE}}
")

  # Define macro for getting PETSc variables from Makefile
  macro(PETSC_GET_VARIABLE var name)
    set(${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE petsc_return)
  endmacro()

  # Call macro to get the PETSc variables
  petsc_get_variable(PETSC_INCLUDE PETSC_INCLUDE)          # 3.1
  petsc_get_variable(PETSC_CC_INCLUDES PETSC_CC_INCLUDES)  # dev
  set(PETSC_INCLUDE ${PETSC_INCLUDE} ${PETSC_CC_INCLUDES})
  petsc_get_variable(PETSC_LIB_BASIC PETSC_LIB_BASIC)
  petsc_get_variable(PETSC_LIB_DIR PETSC_LIB_DIR)
  set(PETSC_LIB "-L${PETSC_LIB_DIR} ${PETSC_LIB_BASIC}")

  # Extract include paths and libraries from compile command line
  include(ResolveCompilerPaths)
  resolve_includes(PETSC_INCLUDE_DIRS "${PETSC_INCLUDE}")
  resolve_libraries(PETSC_LIBRARIES "${PETSC_LIB}")

  # Add some extra libraries on OSX
  if (APPLE)

    # CMake will have troubel finding the gfortan libraries if compiling
    # with clang (the libs may be required by 3rd party Fortran libraries)
    find_program(GFORTRAN_EXECUTABLE gfortran)
    if (GFORTRAN_EXECUTABLE)
      execute_process(COMMAND ${GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.dylib
      OUTPUT_VARIABLE GFORTRAN_LIBRARY
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    find_package(X11)
    list(APPEND PETSC_INCLUDE_DIRS ${X11_X11_INCLUDE_PATH})

  endif()

  # Remove temporary Makefile
  file(REMOVE ${petsc_config_makefile})

  # Add variables to CMake cache and mark as advanced
  set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS} CACHE STRING "PETSc include paths." FORCE)
  set(PETSC_LIBRARIES ${PETSC_LIBRARIES} CACHE STRING "PETSc libraries." FORCE)
  mark_as_advanced(PETSC_INCLUDE_DIRS PETSC_LIBRARIES)

endif()

# Build PETSc test program
if (FOUND_PETSC_CONF)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${PETSC_LIBRARIES})

  # Add MPI variables if MPI has been found
  if (MPI_C_FOUND)
    set(CMAKE_REQUIRED_INCLUDES  ${CMAKE_REQUIRED_INCLUDES} ${MPI_C_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS     "${CMAKE_REQUIRED_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  endif()

  # Check PETSc version
  set(PETSC_CONFIG_TEST_VERSION_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/petsc_config_test_version.cpp")
  file(WRITE ${PETSC_CONFIG_TEST_VERSION_CPP} "
#include <iostream>
#include \"petscversion.h\"

int main() {
  std::cout << PETSC_VERSION_MAJOR << \".\"
	    << PETSC_VERSION_MINOR << \".\"
	    << PETSC_VERSION_SUBMINOR;
  return 0;
}
")

  try_run(
    PETSC_CONFIG_TEST_VERSION_EXITCODE
    PETSC_CONFIG_TEST_VERSION_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PETSC_CONFIG_TEST_VERSION_CPP}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE OUTPUT
    )

  if (PETSC_CONFIG_TEST_VERSION_EXITCODE EQUAL 0)
    set(PETSC_VERSION ${OUTPUT} CACHE TYPE STRING)
    string(REPLACE "." ";" PETSC_VERSION_LIST ${PETSC_VERSION})
    list(GET PETSC_VERSION_LIST 0 PETSC_VERSION_MAJOR)
    list(GET PETSC_VERSION_LIST 1 PETSC_VERSION_MINOR)
    list(GET PETSC_VERSION_LIST 2 PETSC_VERSION_SUBMINOR)
    mark_as_advanced(PETSC_VERSION)
    mark_as_advanced(PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR)
  endif()

  if (PETSc_FIND_VERSION)
    # Check if version found is >= required version
    if (NOT "${PETSC_VERSION}" VERSION_LESS "${PETSc_FIND_VERSION}")
      set(PETSC_VERSION_OK TRUE)
    endif()
  else()
    # No specific version requested
    set(PETSC_VERSION_OK TRUE)
  endif()
  mark_as_advanced(PETSC_VERSION_OK)

else()
  set(PETSC_VERSION "UNKNOWN")
  set(PETSC_VERSION_OK FALSE)
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  "PETSc could not be found. Be sure to set PETSC_DIR and PETSC_ARCH."
  PETSC_LIBRARIES PETSC_DIR PETSC_INCLUDE_DIRS 
  PETSC_VERSION PETSC_VERSION_OK)
