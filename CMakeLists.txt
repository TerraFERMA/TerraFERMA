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

project(terraferma)

set(BUCKETTOOLS_DIR "${PROJECT_SOURCE_DIR}/buckettools")

add_subdirectory(${BUCKETTOOLS_DIR} "buckettools")

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cpp DESTINATION share/terraferma)

file(
    WRITE ${PROJECT_BINARY_DIR}/configmodule
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

install(FILES ${PROJECT_BINARY_DIR}/configmodule DESTINATION share/terraferma)

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

set(CMAKE_MODULE_PATH "${BUCKETTOOLS_DIR}/cmake/modules")
find_package(PETSc)
if(NOT ${PETSc_FOUND})
  message(FATAL_ERROR "Could not find PETSc.")
endif(NOT ${PETSc_FOUND})

set(EXCLUDE_TAGS "not:branch:${GIT_BRANCH} not:petsc:${PETSC_VERSION} not:petsc:${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}")

file(
    WRITE ${PROJECT_BINARY_DIR}/copy_test_input.cmake
"set( LIST_COMMAND tfsimulationharness --list-input -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${LIST_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process (
   COMMAND \${LIST_COMMAND}
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests
   RESULT_VARIABLE RETCODE
   OUTPUT_VARIABLE TEST_INPUT_STRING
   ERROR_VARIABLE TEST_INPUT_ERR
   OUTPUT_STRIP_TRAILING_WHITESPACE
   ERROR_STRIP_TRAILING_WHITESPACE
   )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()

string(REPLACE \"\\n\" \";\" TEST_INPUT_LIST \${TEST_INPUT_STRING})

MACRO(COPY_IF_DIFFERENT source dest)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different \${source} \${dest}
    )
ENDMACRO()

foreach(INPUTFILE \${TEST_INPUT_LIST})
  COPY_IF_DIFFERENT(\"${PROJECT_SOURCE_DIR}/tests/\${INPUTFILE}\" \"${PROJECT_BINARY_DIR}/tests/\${INPUTFILE}\")
endforeach()"
    )

add_custom_target(
  copy_test_input
  COMMAND ${CMAKE_COMMAND} -P copy_test_input.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/clean_test.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( CLEAN_COMMAND tfsimulationharness --clean -n \$ENV{THREADS} -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${CLEAN_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process(
    COMMAND \${CLEAN_COMMAND} 
    RESULT_VARIABLE RETCODE
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()
"
    )

add_custom_target(
  clean_test
  COMMAND ${CMAKE_COMMAND} -P clean_test.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/build_test.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( BUILD_COMMAND tfsimulationharness --build -l short medium long -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${BUILD_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process(
    COMMAND \${BUILD_COMMAND}
    RESULT_VARIABLE RETCODE
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()
"
    )

add_custom_target(
  build_test
  COMMAND ${CMAKE_COMMAND} -P build_test.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/test.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( TEST_COMMAND tfsimulationharness --test -l short -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${TEST_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process(
    COMMAND \${TEST_COMMAND}
    RESULT_VARIABLE RETCODE
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()
"
    )

add_custom_target(
  test
  COMMAND ${CMAKE_COMMAND} -P test.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/mediumtest.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( TEST_COMMAND tfsimulationharness --test -l medium -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${TEST_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process(
    COMMAND \${TEST_COMMAND}
    RESULT_VARIABLE RETCODE
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()
"
    )

add_custom_target(
  mediumtest
  COMMAND ${CMAKE_COMMAND} -P mediumtest.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/longtest.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( TEST_COMMAND tfsimulationharness --test -l long -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo "Running:" \${TEST_COMMAND}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
     )

execute_process(
    COMMAND \${TEST_COMMAND}
    RESULT_VARIABLE RETCODE
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()
"
    )

add_custom_target(
  longtest
  COMMAND ${CMAKE_COMMAND} -P longtest.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

# install tutorials under share
install(DIRECTORY ${PROJECT_SOURCE_DIR}/tutorials DESTINATION share/terraferma)

