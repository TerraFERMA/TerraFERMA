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

find_path(SPUD_INCLUDE_DIRS
  NAMES spud
  HINTS ${SPUD_DIR}/include $ENV{SPUD_DIR}/include
  DOC "Directory where the Spud header file is located"
  )
mark_as_advanced(SPUD_INCLUDE_DIRS)

find_library(SPUD_LIBRARIES
  NAMES spud
  HINTS ${SPUD_DIR}/lib $ENV{SPUD_DIR}/lib
  DOC "The Spud library"
  )
mark_as_advanced(SPUD_LIBRARIES)

set(${SPUD_LIBRARIES} "${SPUD_LIBRARIES}")

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Spud
  "Spud could not be found. Be sure to set SPUD_DIR."
  SPUD_LIBRARIES SPUD_INCLUDE_DIRS)
