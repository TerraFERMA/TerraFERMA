# - Try to find BucketTools
# Once done (and if successful) this will define

#  BUCKETTOOLS_FOUND        - system has BucketTools
#  BUCKETTOOLS_INCLUDE_DIRS - include directories for BucketTools
#  BUCKETTOOLS_LIBRARIES    - libraries for BuckeTools
#  BUCKETTOOLS_SHARED_DIRS  - location of shared material
#  BUCKETTOOLS_BIN          - location of executables
#
# Setting these changes the behavior of the search
#
#  BUCKETTOOLS_DIR - directory in which BucketTools resides

message(STATUS "Checking for package 'BucketTools'")

find_path(BUCKETTOOLS_INCLUDE_DIRS
  NAMES BucketTools.h
  HINTS ${BUCKETTOOLS_DIR}/include $ENV{BUCKETTOOLS_DIR}/include
  DOC "Directory where the BucketTools header file is located"
  )
mark_as_advanced(BUCKETTOOLS_INCLUDE_DIRS)

find_path(BUCKETTOOLS_SHARED_DIRS
  NAMES buckettools
  HINTS ${BUCKETTOOLS_DIR}/share $ENV{BUCKETTOOLS_DIR}/share
  DOC "Directory where the BucketTools shared material is located"
  )
mark_as_advanced(BUCKETTOOLS_SHARED_DIRS)

find_path(BUCKETTOOLS_BIN
  NAMES systemwrappers_from_options
  HINTS ${BUCKETTOOLS_DIR}/bin $ENV{BUCKETTOOLS_DIR}/bin
  DOC "Directory where the BucketTools executables are located"
  )
mark_as_advanced(BUCKETTOOLS_BIN)

find_library(BUCKETTOOLS_LIBRARIES
  NAMES buckettools_cpp
  HINTS ${BUCKETTOOLS_DIR}/lib $ENV{BUCKETTOOLS_DIR}/lib
  DOC "The BucketTools cpp library"
  )
mark_as_advanced(BUCKETTOOLS_LIBRARIES)

set(${BUCKETTOOLS_LIBRARIES} "${BUCKETTOOLS_LIBRARIES}")

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BucketTools
  "BucketTools could not be found. Be sure to set BUCKETTOOLS_DIR."
  BUCKETTOOLS_LIBRARIES BUCKETTOOLS_INCLUDE_DIRS BUCKETTOOLS_SHARED_DIRS BUCKETTOOLS_BIN)

