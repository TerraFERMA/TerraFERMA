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

set(EXECUTABLE "alt_courant_1d")

set(OPTIONSFILENAME "courant.tfml")

set(PREFIX ${EXECUTABLE})

project($EXECUTABLE)

set(OPTIONSFILE "${PROJECT_SOURCE_DIR}/${OPTIONSFILENAME}")

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

add_executable(${EXECUTABLE} "${PROJECT_SOURCE_DIR}/../../cpp/main.cpp")

set (EXTRA_LIBS ${EXTRA_LIBS} ${BUCKETTOOLS_LIBRARIES} ${PREFIX}buckettools_ufc)

target_link_libraries(${EXECUTABLE} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${EXTRA_LIBS} ${SPUD_LIBRARIES} ${PYTHON_LIBRARIES})

add_custom_target(
  ${PREFIX}_run
  COMMAND ${EXECUTABLE} -v${LOGLEVEL} -l ${OPTIONSFILE} 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${EXECUTABLE}
  )

add_custom_target(
  ${PREFIX}_test
  COMMAND py.test -vs ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )
add_dependencies(${PREFIX}_test ${PREFIX}_run)
