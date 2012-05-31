
add_subdirectory(cpp)
add_subdirectory(ufc)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

if (NOT DEFINED BUCKETTOOLS_PYTHON_INSTALL_DIR)
  # Get path for pure Python modules
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE BUCKETTOOLS_PYTHON_INSTALL_DIR
    )
  # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
  string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
    BUCKETTOOLS_PYTHON_INSTALL_DIR "${BUCKETTOOLS_PYTHON_INSTALL_DIR}")
  set(BUCKETTOOLS_PYTHON_INSTALL_DIR ${BUCKETTOOLS_PYTHON_INSTALL_DIR}
    CACHE PATH "Python module installation directory.")
endif()

if (BUCKETTOOLS_PYTHON_INSTALL_DIR)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/buckettools DESTINATION ${BUCKETTOOLS_PYTHON_INSTALL_DIR})
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/schemas DESTINATION share/buckettools)

# Write a schemata file for bmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/bml
    "BucketTools Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/buckettools.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/bml DESTINATION share/buckettools/diamond/schemata)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake DESTINATION share/buckettools)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/ufc   DESTINATION share/buckettools)

