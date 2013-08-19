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

set(PROJECT_NAME advection)
project(${PROJECT_NAME})

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES})
