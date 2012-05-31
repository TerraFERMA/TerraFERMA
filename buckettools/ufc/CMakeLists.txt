# The project name for buckettools (sets up environment variables for binary and source directories)
project(buckettools_ufc)

# add a command that runs the ufc generation script if the options file
# or the (main file of the) script has changed
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp ${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp ${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp
    COMMAND ${BUCKETTOOLS_BIN}/systemwrappers_from_options
    ARGS ${OPTIONSFILE}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${OPTIONSFILE}
    DEPENDS ${BUCKETTOOLS_BIN}/systemwrappers_from_options
    )
 
# create a custom target so we can make the project depend on this
add_custom_target(
   ${PREFIX}genwrappers
   DEPENDS ${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp ${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp ${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp
   )

# include the buckettools include directory and the binary directory (for the systemswrapper)
include_directories(${BUCKETTOOLS_INCLUDE_DIRS} ${PROJECT_BINARY_DIR})
# include python
include_directories(${PYTHON_INCLUDE_PATH})
# generate a library - really the main point of this whole process
add_library(${PREFIX}buckettools_ufc ${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp ${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp ${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp)
# tell cmake that this file doesn't exist until build time
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemFunctionalsWrapper.cpp PROPERTIES GENERATED 1)
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemSolversWrapper.cpp PROPERTIES GENERATED 1)
set_source_files_properties(${PROJECT_BINARY_DIR}/SystemExpressionsWrapper.cpp PROPERTIES GENERATED 1)
# the project depends on this target
add_dependencies(${PREFIX}buckettools_ufc ${PREFIX}genwrappers)
# link to other libraries
target_link_libraries(${PREFIX}buckettools_ufc ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES})

