project(buckettools)

file(
    WRITE ${buckettools_BINARY_DIR}/confdefs.h.in
    "#define __HG_ID__ \"@HG_ID@\"\n"
    )

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

add_custom_target(
    confdefs
    ${CMAKE_COMMAND} -D SRC="${buckettools_BINARY_DIR}/confdefs.h.in"
                     -D DST="${buckettools_BINARY_DIR}/confdefs.h"
                     -D DIR="${PROJECT_SOURCE_DIR}"
                     -P "${buckettools_BINARY_DIR}/confdefs.cmake"
    )

add_subdirectory("ufl")

include_directories("${buckettools_SOURCE_DIR}/include" ${buckettools_BINARY_DIR} ${buckettools_BINARY_DIR}/ufl)
include_directories(${PYTHON_INCLUDE_PATH})
add_library(buckettools cpp/Bucket.cpp cpp/SpudBucket.cpp cpp/System.cpp cpp/SpudSystem.cpp)
set_source_files_properties(${buckettools_BINARY_DIR}/ufl/SystemsWrapper.h PROPERTIES GENERATED 1)
add_dependencies(buckettools generate_ufl)
set_source_files_properties(confdefs.h PROPERTIES GENERATED 1)
add_dependencies(buckettools confdefs)
target_link_libraries(buckettools ${PYTHON_LIBRARIES})


