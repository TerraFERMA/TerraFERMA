#file(
#    WRITE ${BucketTools_BINARY_DIR}/confdefs.h.in
#    "#define __HG_ID__ \"@HG_ID@\"\n"
#    )
#
#file(
#    WRITE ${BucketTools_BINARY_DIR}/confdefs.cmake
#    "find_program(HAVE_HG hg)
#    if(HAVE_HG)
#        execute_process (
#           COMMAND hg id -i 
#           WORKING_DIRECTORY \${DIR}
#           OUTPUT_VARIABLE HG_ID
#           ERROR_VARIABLE HG_ERR
#           OUTPUT_STRIP_TRAILING_WHITESPACE
#           ERROR_STRIP_TRAILING_WHITESPACE
#           )
#        if(HG_ERR)
#           set(HG_ID "hg_id_not_found")
#        endif(HG_ERR)
#    else(HAVE_HG)
#        set(HG_ID "hg_not_found")
#    endif(HAVE_HG)
#    configure_file(\${SRC} \${DST} @ONLY)
#    "
#    )
#
#add_custom_target(
#    confdefs
#    ${CMAKE_COMMAND} -D SRC="${BucketTools_BINARY_DIR}/confdefs.h.in"
#                     -D DST="${BucketTools_BINARY_DIR}/confdefs.h"
#                     -D DIR="${PROJECT_SOURCE_DIR}"
#                     -P "${BucketTools_BINARY_DIR}/confdefs.cmake"
#    )
#
#include_directories("${BucketTools_SOURCE_DIR}/include" ${BucketTools_BINARY_DIR})
#include_directories(${PYTHON_INCLUDE_PATH})
#add_library(BucketTools Bucket.cpp SpudBucket.cpp System.cpp SpudSystem.cpp)
##add_library(BucketTools Bucket.cpp SystemBucket.cpp SpudBucket.cpp
##                        StatFile.cpp DiagnosticsFile.cpp DetectorsFile.cpp
##                        PythonExpression.cpp InitialConditionExpression.cpp PythonInstance.cpp 
##                        GenericDetectors.cpp PointDetectors.cpp PythonDetectors.cpp)
#add_dependencies(BucketTools confdefs)
##target_link_libraries(BucketTools BucketUFCFiles ${PYTHON_LIBRARIES})
#target_link_libraries(BucketTools ${PYTHON_LIBRARIES})
#

