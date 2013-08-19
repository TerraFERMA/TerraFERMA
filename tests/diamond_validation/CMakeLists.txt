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

set(EXECUTABLE "diamond_validation")

set(PREFIX ${EXECUTABLE})

project($EXECUTABLE)

add_custom_target(
  ${PREFIX}_run
  COMMAND python ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )

add_custom_target(
  ${PREFIX}_test
  COMMAND py.test -vs ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )

