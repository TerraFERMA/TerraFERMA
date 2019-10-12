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

file(
    WRITE ${PROJECT_BINARY_DIR}/copy_test_input.cmake
"set( LIST_COMMAND tfsimulationharness --list-input -r -- *.shml )

execute_process(
     COMMAND echo Running: \${LIST_COMMAND}
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
    WRITE ${PROJECT_BINARY_DIR}/clean_tests.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( CLEAN_COMMAND tfsimulationharness --clean -n \$ENV{THREADS} -r -- *.shml )

execute_process(
     COMMAND echo Running: \${CLEAN_COMMAND}
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
  clean_tests
  COMMAND ${CMAKE_COMMAND} -P clean_tests.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/build_tests.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

set( BUILD_COMMAND tfsimulationharness --build -l short medium long -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo Running: \${BUILD_COMMAND}
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
  build_tests
  COMMAND ${CMAKE_COMMAND} -P build_tests.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/run_shorttests.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

if (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --run -l short --mpi-options ' --oversubscribe' -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
else (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --run -l short -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
endif (ENV{OVERSUBSCRIBE})

execute_process(
     COMMAND echo Running: \${TEST_COMMAND}
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

if (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --run-from-checkpoint -l short --mpi-options ' --oversubscribe' -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -t run_from_checkpoint -r -- *.shml )
else (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --run-from-checkpoint -l short -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -t run_from_checkpoint -r -- *.shml )
endif (ENV{OVERSUBSCRIBE})

execute_process(
     COMMAND echo Running: \${TEST_COMMAND}
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

set( TEST_COMMAND tfsimulationharness --just-test -l short -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )

execute_process(
     COMMAND echo Running: \${TEST_COMMAND}
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
  run_shorttests
  COMMAND ${CMAKE_COMMAND} -P run_shorttests.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/run_mediumtests.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

if (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --test -l medium --mpi-options ' --oversubscribe' -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
else (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --test -l medium -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
endif (ENV{OVERSUBSCRIBE})

execute_process(
     COMMAND echo Running: \${TEST_COMMAND}
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
  run_mediumtests
  COMMAND ${CMAKE_COMMAND} -P run_mediumtests.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

file(
    WRITE ${PROJECT_BINARY_DIR}/run_longtests.cmake
"if (NOT DEFINED ENV{THREADS})
  set( ENV{THREADS} 1 )
endif()

if (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --test -l long --mpi-options ' --oversubscribe' -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
else (ENV{OVERSUBSCRIBE})
  set( TEST_COMMAND tfsimulationharness --test -l long -n \$ENV{THREADS} -e ${EXCLUDE_TAGS} -r -- *.shml )
endif (ENV{OVERSUBSCRIBE})

execute_process(
     COMMAND echo Running: ${TEST_COMMAND}
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
  run_longtests
  COMMAND ${CMAKE_COMMAND} -P run_longtests.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS copy_test_input
  )

# install tests under share
install(DIRECTORY ${PROJECT_SOURCE_DIR}/tests DESTINATION share/terraferma)

