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
project(buckettools_cpp)

# Write a template header file to disk at configuration time
# (this holds the template for the mercurial id information that will
#  be replaced at build time by a suitable variable - using the script
#  below)
file(
    WRITE ${PROJECT_BINARY_DIR}/builddefs.h.in
    "#define __GIT_SHA__ @GIT_SHA@\n"
    )

# Write a cmake script to disk at configuration time
# (this will run git - at compile time - and record the id info
#  to an environment variable which can be used to write a
#  header file for the build)
file(
    WRITE ${PROJECT_BINARY_DIR}/builddefs.cmake
    "find_program(HAVE_GIT git)
    if(HAVE_GIT)
        execute_process (
           COMMAND git log --pretty=format:\"%H %ad\" -1 HEAD
           WORKING_DIRECTORY \${DIR}
           OUTPUT_VARIABLE GIT_SHA
           ERROR_VARIABLE GIT_ERR
           OUTPUT_STRIP_TRAILING_WHITESPACE
           ERROR_STRIP_TRAILING_WHITESPACE
           )
        if(GIT_ERR)
           set(GIT_SHA git_log_not_found)
        endif(GIT_ERR)
    else(HAVE_GIT)
        set(GIT_SHA git_not_found)
    endif(HAVE_GIT)
    configure_file(\${SRC} \${DST} @ONLY)
    "
    )

# add a target that runs the above cmake script and creates the header
# file containing the git log information
add_custom_target(
    builddefs
    ${CMAKE_COMMAND} -D SRC="${PROJECT_BINARY_DIR}/builddefs.h.in"
                     -D DST="${PROJECT_BINARY_DIR}/builddefs.h"
                     -D DIR="${PROJECT_SOURCE_DIR}"
                     -P "${PROJECT_BINARY_DIR}/builddefs.cmake"
    )


# include the buckettools include directory and the binary directory (for builddefs and the systemswrapper)
include_directories("${PROJECT_SOURCE_DIR}/../include" ${PROJECT_BINARY_DIR})
# include python
include_directories(${PYTHON_INCLUDE_PATH})
# generate a library - really the main point of this whole process
add_library(buckettools_cpp SHARED
                            Usage.cpp SignalHandler.cpp SigIntEventHandler.cpp
                            Bucket.cpp SpudBucket.cpp SystemBucket.cpp SpudSystemBucket.cpp
                            FunctionBucket.cpp SpudFunctionBucket.cpp
                            SolverBucket.cpp SpudSolverBucket.cpp 
                            SpudBase.cpp MPIBase.cpp PythonExpression.cpp PythonInstance.cpp
                            InitialConditionExpression.cpp RegionsExpression.cpp SemiLagrangianExpression.cpp
                            GenericDetectors.cpp PointDetectors.cpp PythonDetectors.cpp
                            DiagnosticsFile.cpp StatisticsFile.cpp SteadyStateFile.cpp
                            DetectorsFile.cpp ConvergenceFile.cpp KSPConvergenceFile.cpp
                            BucketPETScBase.cpp BucketDolfinBase.cpp DolfinPETScBase.cpp
                            ReferencePoints.cpp)
# tell cmake that this file doesn't exist until build time
set_source_files_properties(builddefs.h PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(buckettools_cpp builddefs)
# link to other libraries
target_link_libraries(buckettools_cpp ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${SPUD_LIBRARIES} ${PYTHON_LIBRARIES})
get_target_property(BUCKETTOOLS_LIBRARY_LOCATION buckettools_cpp LOCATION)
get_filename_component(BUCKETTOOLS_LIBRARY_FILENAME ${BUCKETTOOLS_LIBRARY_LOCATION} NAME)
set(BUCKETTOOLS_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/${BUCKETTOOLS_LIBRARY_FILENAME}" PARENT_SCOPE)

install(TARGETS buckettools_cpp DESTINATION lib)


