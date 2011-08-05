# The project name for buckettools (sets up environment variables for binary and source directories)
project(buckettools_ufc)

# add a command that runs the ufc generation script if the options file
# or the (main file of the) script has changed
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/SystemsWrapper.cpp
    COMMAND python
    ARGS ${PROJECT_SOURCE_DIR}/../scripts/ufc_from_options.py ${OPTIONSFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${OPTIONSFILE}
    DEPENDS ${PROJECT_SOURCE_DIR}/../scripts/ufc_from_options.py
    )
 
# create a custom target so we can make the project depend on this
add_custom_target(
   ${PREFIX}genufc
   DEPENDS ${PROJECT_BINARY_DIR}/SystemsWrapper.cpp
   )

# include the buckettools include directory and the binary directory (for the systemswrapper)
include_directories("${PROJECT_SOURCE_DIR}/../include" ${PROJECT_BINARY_DIR})
# include python
include_directories(${PYTHON_INCLUDE_PATH})
# generate a library - really the main point of this whole process
add_library(${PREFIX}buckettools_ufc ${PROJECT_BINARY_DIR}/SystemsWrapper.cpp )
# tell cmake that this file doesn't exist until build time
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemsWrapper.cpp PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(${PREFIX}buckettools_ufc ${PREFIX}genufc)

