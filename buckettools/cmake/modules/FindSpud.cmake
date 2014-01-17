# - Try to find Spud
# Once done this will define

#  SPUD_FOUND        - system has Spud
#  SPUD_INCLUDE_DIRS - include directories for Spud
#  SPUD_LIBRARIES    - libraries for Spud
#
# Setting these changes the behavior of the search
#
#  SPUD_DIR - directory in which Spud resides

message(STATUS "Checking for package 'Spud'")

find_path(SPUD_INCLUDE_DIRS spud.h
  PATHS ${SPUD_DIR}/include $ENV{SPUD_DIR}/include
  DOC "Directory where the Spud header file is located"
  NO_DEFAULT_PATH
  )
find_path(SPUD_INCLUDE_DIRS spud.h)
mark_as_advanced(SPUD_INCLUDE_DIRS)

# Report result of search for SPUD_INCLUDE_DIRS
message(STATUS "SPUD_INCLUDE_DIRS is ${SPUD_INCLUDE_DIRS}")

find_library(SPUD_LIBRARIES spud
  PATHS ${SPUD_DIR}/lib $ENV{SPUD_DIR}/lib
  DOC "The Spud library"
  NO_DEFAULT_PATH
  )
find_library(SPUD_LIBRARIES spud)
mark_as_advanced(SPUD_LIBRARIES)

# Report result of search for SPUD_LIBRARIES
message(STATUS "SPUD_LIBRARIES is ${SPUD_LIBRARIES}")

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Spud
  "Spud could not be found. Be sure to set SPUD_DIR."
  SPUD_LIBRARIES SPUD_INCLUDE_DIRS)
