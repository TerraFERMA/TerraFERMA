
add_custom_command(
    OUTPUT ${BucketTools_BINARY_DIR}/ufl/SystemsWrapper.h
    COMMAND python
    ARGS ${BucketTools_SOURCE_DIR}/scripts/ufl_from_options.py ${OPTIONSFILE}
    WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
    DEPENDS ${OPTIONSFILE}
    DEPENDS ${BucketTools_SOURCE_DIR}/scripts/ufl_from_options.py
    )

add_custom_target(
    generate_ufl 
    DEPENDS ${BucketTools_BINARY_DIR}/ufl/SystemsWrapper.h
    )

