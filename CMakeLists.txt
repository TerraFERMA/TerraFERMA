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

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cpp DESTINATION share)

file(
    WRITE ${PROJECT_BINARY_DIR}/basebucket-config
"#%Module1.0
proc ModulesHelp { } {
        global dotversion

        puts stderr \"\tSets environment for Basebucket.\"
}
module-whatis \"Sets environment for Basebucket and Buckettools.\"

setenv       BASEBUCKET_CMAKE_PATH ${CMAKE_INSTALL_PREFIX}/share/cpp
setenv       BUCKETTOOLS_DIR       ${CMAKE_INSTALL_PREFIX}
prepend-path CMAKE_MODULE_PATH     ${CMAKE_INSTALL_PREFIX}/share/buckettools/cmake/modules
prepend-path PYTHONPATH            ${CMAKE_INSTALL_PREFIX}/${BUCKETTOOLS_PYTHON_INSTALL_DIR}
prepend-path DIAMOND_CONFIG_PATH   ${CMAKE_INSTALL_PREFIX}/share/buckettools/diamond
prepend-path PATH                  ${CMAKE_INSTALL_PREFIX}/bin"
    )

install(FILES ${PROJECT_BINARY_DIR}/basebucket-config DESTINATION share)




