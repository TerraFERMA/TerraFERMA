
add_custom_command(
    OUTPUT ${buckettools_BINARY_DIR}/ufl/SystemsWrapper.h
    COMMAND python
    ARGS ${buckettools_SOURCE_DIR}/scripts/ufl_from_options.py ${OPTIONSFILE}
    WORKING_DIRECTORY ${buckettools_BINARY_DIR}/ufl
    DEPENDS ${OPTIONSFILE}
    DEPENDS ${buckettools_SOURCE_DIR}/scripts/ufl_from_options.py
    )

add_custom_target(
    generate_ufl 
    DEPENDS ${buckettools_BINARY_DIR}/ufl/SystemsWrapper.h
    )

