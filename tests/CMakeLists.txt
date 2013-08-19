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

project(tests)

if (NOT LOGLEVEL)
  set(LOGLEVEL "INFO" CACHE STRING
       "The log level that should be used in run targets, options are: DBG, TRACE, PROGRESS, INFO, WARNING, ERROR, CRITICAL." FORCE)
endif()

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN)

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

find_package(PythonLibs REQUIRED)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH} "${BUCKETTOOLS_DIR}/cmake/modules")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_BINARY_DIR}/buckettools)

find_package(BucketTools REQUIRED)

if (NOT BUCKETTOOLS_FOUND)
  message(FATAL_ERROR "Could not find BucketTools.")
endif()

include_directories(${BUCKETTOOLS_INCLUDE_DIRS})

include_directories(${PYTHON_INCLUDE_PATH})

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

add_subdirectory(diamond_validation)
add_subdirectory(courant_1d)
add_subdirectory(alt_courant_1d)
add_subdirectory(adaptive_dt_1d)
add_subdirectory(alt_adaptive_dt_1d)
add_subdirectory(poisson)
add_subdirectory(neumann_poisson)
add_subdirectory(projection_mixedelement)
add_subdirectory(projection_timedependent)
add_subdirectory(projection_timedependent_checkpoint)
add_subdirectory(projection_cpp)
add_subdirectory(rbconvection_steadystate_split)
add_subdirectory(rbconvection_steadystate_split_schur)
add_subdirectory(rbconvection_steadystate)
add_subdirectory(rbconvection_steadystate_schur)
add_subdirectory(rbconvection1a)
add_subdirectory(rbconvection1a_refpt)
add_subdirectory(rbconvection2a)
add_subdirectory(region_ids)
add_subdirectory(subductionbenchmark1a)
add_subdirectory(subductionbenchmark1b)
add_subdirectory(subductionbenchmark1c)
add_subdirectory(subductionbenchmark1c_split)
add_subdirectory(subductionbenchmark2a)
add_subdirectory(subductionbenchmark2a_split)
add_subdirectory(subductionbenchmark2b)
add_subdirectory(dg_advection_1d_vi)
add_subdirectory(viscosity_contrast)

add_custom_target(build_test)
add_dependencies(build_test courant_1d alt_courant_1d adaptive_dt_1d alt_adaptive_dt_1d
projection_mixedelement projection_timedependent projection_timedependent_checkpoint projection_cpp
rbconvection_steadystate_split rbconvection_steadystate_split_schur rbconvection_steadystate 
rbconvection_steadystate_schur rbconvection1a rbconvection1a_refpt region_ids 
subductionbenchmark1a subductionbenchmark1b subductionbenchmark1c subductionbenchmark1c_split
subductionbenchmark2a subductionbenchmark2a_split subductionbenchmark2b neumann_poisson 
dg_advection_1d_vi viscosity_contrast)

add_custom_target(run)
add_dependencies(run poisson_run)

add_custom_target(run3p1)
add_dependencies(run3p1 run)

add_custom_target(run3p2)
add_dependencies(run3p2 run)

add_custom_target(run3p3)
add_dependencies(run3p3 run)

add_custom_target(run_test3p1)
add_dependencies(run_test3p1 diamond_validation_test courant_1d_test alt_courant_1d adaptive_dt_1d_test alt_adaptive_dt_1d_test
projection_mixedelement_test projection_timedependent_test projection_timedependent_checkpoint_test projection_cpp_test
rbconvection_steadystate_split_test rbconvection_steadystate_test 
rbconvection1a_test rbconvection1a_refpt_test region_ids_test 
subductionbenchmark1a_test subductionbenchmark1b_test subductionbenchmark1c_test subductionbenchmark1c_split_test
subductionbenchmark2a_test subductionbenchmark2a_split_test subductionbenchmark2b_test neumann_poisson_test)

add_custom_target(run_test3p2)
add_dependencies(run_test3p2 run_test3p1 dg_advection_1d_vi_test rbconvection_steadystate_split_schur_test
rbconvection_steadystate_schur_test)

add_custom_target(run_test3p3)
add_dependencies(run_test3p3 run_test3p1 dg_advection_1d_vi_test rbconvection_steadystate_split_schur_test
rbconvection_steadystate_schur_test viscosity_contrast_test)

add_custom_target(run_mediumtest)
add_dependencies(run_mediumtest rbconvection2a_test)

