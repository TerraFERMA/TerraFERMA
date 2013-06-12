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
           set(GIT_SHA "git_log_not_found")
        endif(GIT_ERR)
    else(HAVE_GIT)
        set(GIT_SHA "git_not_found")
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
                            SpudBase.cpp PythonExpression.cpp PythonInstance.cpp
                            InitialConditionExpression.cpp RegionsExpression.cpp SemiLagrangianExpression.cpp
                            GenericDetectors.cpp PointDetectors.cpp PythonDetectors.cpp
                            DiagnosticsFile.cpp StatisticsFile.cpp SteadyStateFile.cpp
                            DetectorsFile.cpp ConvergenceFile.cpp KSPConvergenceFile.cpp
                            BucketPETScBase.cpp BucketDolfinBase.cpp
                            ReferencePoints.cpp)
# tell cmake that this file doesn't exist until build time
set_source_files_properties(builddefs.h PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(buckettools_cpp builddefs)
# link to other libraries
target_link_libraries(buckettools_cpp ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${SPUD_LIBRARIES} ${PYTHON_LIBRARIES})

install(TARGETS buckettools_cpp DESTINATION lib)


