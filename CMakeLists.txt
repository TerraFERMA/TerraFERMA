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

if (NOT BUCKETTOOLS_DIR)
  set(BUCKETTOOLS_DIR "${PROJECT_SOURCE_DIR}/buckettools" CACHE PATH
       "Location of buckettools." FORCE)
endif()

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
setenv       BUCKETTOOLS_DIR       ${CMAKE_INSTALL_PREFIX}
prepend-path CMAKE_MODULE_PATH     ${CMAKE_INSTALL_PREFIX}/share/buckettools/cmake/modules
prepend-path PYTHONPATH            ${CMAKE_INSTALL_PREFIX}/${BUCKETTOOLS_PYTHON_INSTALL_DIR}
prepend-path DIAMOND_CONFIG_PATH   ${CMAKE_INSTALL_PREFIX}/share/buckettools/diamond
prepend-path PATH                  ${CMAKE_INSTALL_PREFIX}/bin

# Special Environment variables for resolving Dynamic Libraries on MacOSX
prepend-path DYLD_LIBRARY_PATH     ${CMAKE_INSTALL_PREFIX}/lib
"

    )

install(FILES ${PROJECT_BINARY_DIR}/configmodule DESTINATION share/terraferma)

option(BUCKET_ENABLE_TESTING "Enable tests." OFF)

if (BUCKET_ENABLE_TESTING)
  add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()


# install tutorials under share
install(DIRECTORY ${PROJECT_SOURCE_DIR}/tutorials DESTINATION share/terraferma)

