# The project name for buckettools (sets up environment variables for binary and source directories)
project(buckettools)

# Write a template header file to disk at configuration time
# (this holds the template for the mercurial id information that will
#  be replaced at build time by a suitable variable - using the script
#  below)
file(
    WRITE ${buckettools_BINARY_DIR}/confdefs.h.in
    "#define __HG_ID__ \"@HG_ID@\"\n"
    )

# Write a cmake script to disk at configuration time
# (this will run hg - at compile time - and record the id info
#  to an environment variable which can be used to write a
#  header file for the build)
file(
    WRITE ${buckettools_BINARY_DIR}/confdefs.cmake
    "find_program(HAVE_HG hg)
    if(HAVE_HG)
        execute_process (
           COMMAND hg id -i 
           WORKING_DIRECTORY \${DIR}
           OUTPUT_VARIABLE HG_ID
           ERROR_VARIABLE HG_ERR
           OUTPUT_STRIP_TRAILING_WHITESPACE
           ERROR_STRIP_TRAILING_WHITESPACE
           )
        if(HG_ERR)
           set(HG_ID "hg_id_not_found")
        endif(HG_ERR)
    else(HAVE_HG)
        set(HG_ID "hg_not_found")
    endif(HAVE_HG)
    configure_file(\${SRC} \${DST} @ONLY)
    "
    )

# add a target that runs the above cmake script and creates the header
# file containing the hg id information
add_custom_target(
    confdefs
    ${CMAKE_COMMAND} -D SRC="${buckettools_BINARY_DIR}/confdefs.h.in"
                     -D DST="${buckettools_BINARY_DIR}/confdefs.h"
                     -D DIR="${PROJECT_SOURCE_DIR}/../"
                     -P "${buckettools_BINARY_DIR}/confdefs.cmake"
    )


# add a command that runs the ufc generation script if the options file
# or the (main file of the) script has changed
add_custom_command(
    OUTPUT ${buckettools_BINARY_DIR}/SystemsWrapper.cpp
    COMMAND python
    ARGS ${buckettools_SOURCE_DIR}/scripts/ufc_from_options.py ${OPTIONSFILE}
    WORKING_DIRECTORY ${buckettools_BINARY_DIR}
    DEPENDS ${OPTIONSFILE}
    DEPENDS ${buckettools_SOURCE_DIR}/scripts/ufc_from_options.py
    )
 
# create a custom target so we can make the project depend on this
add_custom_target(
   genufc
   DEPENDS ${buckettools_BINARY_DIR}/SystemsWrapper.cpp
   )

# include the buckettools include directory and the binary directory (for confdefs and the systemswrapper)
include_directories("${buckettools_SOURCE_DIR}/include" ${buckettools_BINARY_DIR})
# include python
include_directories(${PYTHON_INCLUDE_PATH})
# generate a library - really the main point of this whole process
add_library(buckettools cpp/Bucket.cpp cpp/SpudBucket.cpp cpp/SystemBucket.cpp cpp/SpudSystemBucket.cpp
                        cpp/FunctionBucket.cpp cpp/SpudFunctionBucket.cpp
                        cpp/SolverBucket.cpp cpp/SpudSolverBucket.cpp 
                        cpp/SpudBase.cpp cpp/PythonExpression.cpp cpp/PythonInstance.cpp
                        cpp/InitialConditionExpression.cpp
                        cpp/StatFile.cpp cpp/DiagnosticsFile.cpp
                        cpp/BucketPETScBase.cpp
                        ${buckettools_BINARY_DIR}/SystemsWrapper.cpp )
# tell cmake that this file doesn't exist until build time
set_source_files_properties(${buckettools_BINARY_DIR}/SystemsWrapper.cpp PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(buckettools genufc)
# tell cmake that this file doesn't exist until build time
set_source_files_properties(confdefs.h PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(buckettools confdefs)
# link to python
target_link_libraries(buckettools ${PYTHON_LIBRARIES})


