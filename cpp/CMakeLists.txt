# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(dolfin)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

project(basebucket)

if (NOT EXECUTABLE)
  set(EXECUTABLE "basebucket" CACHE STRING
       "The name of the executable" FORCE)
endif()

if (NOT LOGLEVEL)
  set(LOGLEVEL "INFO" CACHE STRING
       "The log level used in the run target, options are: DBG, TRACE, PROGRESS, INFO, WARNING, ERROR, CRITICAL." FORCE)
endif()

set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})

find_package(BucketTools REQUIRED)

if (NOT BUCKETTOOLS_FOUND)
  message(FATAL_ERROR "Could not find BucketTools.")
endif()

include_directories(${BUCKETTOOLS_INCLUDE_DIRS})

find_package(PythonLibs REQUIRED)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

include_directories(${PYTHON_INCLUDE_PATH})

find_package(PythonInterp REQUIRED)

if(NOT ${PYTHONINTERP_FOUND})
  message(FATAL_ERROR "Could not find python interpreter.")
endif(NOT ${PYTHONINTERP_FOUND})

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

set(OPTIONSFILE CACHE FILEPATH "Choose the input file.")

if(NOT IS_ABSOLUTE ${OPTIONSFILE})
  set(OPTIONSFILE ${PROJECT_BINARY_DIR}/${OPTIONSFILE})
endif(NOT IS_ABSOLUTE ${OPTIONSFILE})

if (IS_DIRECTORY ${OPTIONSFILE})
  message(FATAL_ERROR "OPTIONSFILE not set or set to a directory.")
endif(IS_DIRECTORY ${OPTIONSFILE})

if (NOT EXISTS ${OPTIONSFILE})
  message(FATAL_ERROR "Could not find OPTIONSFILE ${OPTIONSFILE}")
endif(NOT EXISTS ${OPTIONSFILE})

add_subdirectory(${BUCKETTOOLS_SHARED_DIRS}/buckettools/ufc "buckettools_ufc")

add_executable(${EXECUTABLE} ${PROJECT_SOURCE_DIR}/main.cpp)

set (EXTRA_LIBS ${EXTRA_LIBS} ${BUCKETTOOLS_LIBRARIES} buckettools_ufc)

target_link_libraries(${EXECUTABLE} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${EXTRA_LIBS} ${SPUD_LIBRARIES} ${PYTHON_LIBRARIES})

add_custom_target(
  run
  COMMAND ${EXECUTABLE} -v${LOGLEVEL} -l ${OPTIONSFILE} 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${EXECUTABLE}
  )

