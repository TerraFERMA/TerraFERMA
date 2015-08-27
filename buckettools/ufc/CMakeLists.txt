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

# The project name for buckettools (sets up environment variables for binary and source directories)
project(BUCKETTOOLS_UFC)

## write a script that runs the ufc generation script and reruns cmake if necessary
file(
    WRITE ${PROJECT_BINARY_DIR}/genwrappers.cmake
"execute_process(
    COMMAND ${BUCKETTOOLS_BIN}/systemwrappers_from_options -l ${PROJECT_BINARY_DIR}/cpp_filenames.txt.temp -- ${OPTIONSFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    RESULT_VARIABLE RETCODE
    )

if (NOT RETCODE EQUAL 0)
  message(FATAL_ERROR \"Command returned \${RETCODE}\")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files ${PROJECT_BINARY_DIR}/cpp_filenames.txt.temp ${PROJECT_BINARY_DIR}/cpp_filenames.txt
    RESULT_VARIABLE CPP_FILENAMES_CHANGED
    OUTPUT_QUIET
    ERROR_QUIET
    )

if (CPP_FILENAMES_CHANGED)
  execute_process(
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/cpp_filenames.txt.temp ${PROJECT_BINARY_DIR}/cpp_filenames.txt
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )

  execute_process(
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target rebuild_cache
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
endif (CPP_FILENAMES_CHANGED)
"
    )

# create a custom target so we can make the project depend on this
add_custom_target(
   genwrappers
   COMMAND ${CMAKE_COMMAND} -P genwrappers.cmake
   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
   )

if(NOT EXISTS ${PROJECT_BINARY_DIR}/cpp_filenames.txt)
  # perform a dry run of code generation so that a base cpp_filenames can be created at configuration time
  execute_process(
     COMMAND ${BUCKETTOOLS_BIN}/systemwrappers_from_options -d -l ${PROJECT_BINARY_DIR}/cpp_filenames.txt -- ${OPTIONSFILE}
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
     )
endif()

file(STRINGS ${PROJECT_BINARY_DIR}/cpp_filenames.txt CPP_FILENAMES)

# generate a library - really the main point of this whole process
add_library(buckettools_ufc SHARED
            ${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp ${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp ${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp
            ${CPP_FILENAMES})
# tell cmake that this file doesn't exist until build time
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp PROPERTIES GENERATED 1)
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp PROPERTIES GENERATED 1)
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp PROPERTIES GENERATED 1)
set_source_files_properties(${CPP_FILENAMES} PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(buckettools_ufc genwrappers)
# link to other libraries
target_link_libraries(buckettools_ufc ${DOLFIN_LIBRARIES})

