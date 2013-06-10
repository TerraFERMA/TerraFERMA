# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

project(basebucket)

if (NOT BUCKETTOOLS_DIR)
  set(BUCKETTOOLS_DIR "${PROJECT_SOURCE_DIR}/buckettools" CACHE PATH
       "Location of buckettools." FORCE)
endif()

add_subdirectory(${BUCKETTOOLS_DIR} "buckettools")

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cpp DESTINATION share/basebucket)

file(
    WRITE ${PROJECT_BINARY_DIR}/configmodule
"#%Module1.0
proc ModulesHelp { } {
        global dotversion

        puts stderr \"\tSets environment for Basebucket.\"
}
module-whatis \"Sets environment for Basebucket and Buckettools.\"

setenv       BASEBUCKET_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/share/basebucket/cpp
setenv       BUCKETTOOLS_DIR       ${CMAKE_INSTALL_PREFIX}
prepend-path CMAKE_MODULE_PATH     ${CMAKE_INSTALL_PREFIX}/share/buckettools/cmake/modules
prepend-path PYTHONPATH            ${CMAKE_INSTALL_PREFIX}/${BUCKETTOOLS_PYTHON_INSTALL_DIR}
prepend-path DIAMOND_CONFIG_PATH   ${CMAKE_INSTALL_PREFIX}/share/buckettools/diamond
prepend-path PATH                  ${CMAKE_INSTALL_PREFIX}/bin

# Special Environment variables for resolving Dynamic Libraries on MacOSX
prepend-path DYLD_LIBRARY_PATH     ${CMAKE_INSTALL_PREFIX}/lib
"

    )

install(FILES ${PROJECT_BINARY_DIR}/configmodule DESTINATION share/basebucket)

option(BUCKET_ENABLE_TESTING "Enable tests." OFF)

if (BUCKET_ENABLE_TESTING)
  add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()


