set(EXECUTABLE "projection")

set(OPTIONSFILENAME "${EXECUTABLE}.bml")

set(PREFIX ${EXECUTABLE})

project($EXECUTABLE)

set(OPTIONSFILE "${PROJECT_SOURCE_DIR}/${OPTIONSFILENAME}")

if(NOT IS_ABSOLUTE ${OPTIONSFILE})
  set(OPTIONSFILE ${PROJECT_BINARY_DIR}/${OPTIONSFILE})
endif(NOT IS_ABSOLUTE ${OPTIONSFILE})

if (IS_DIRECTORY ${OPTIONSFILE})
  message(FATAL_ERROR "OPTIONSFILE not set or set to a directory.")
endif(IS_DIRECTORY ${OPTIONSFILE})

if (NOT EXISTS ${OPTIONSFILE})
  message(FATAL_ERROR "Could not find OPTIONSFILE ${OPTIONSFILE}")
endif(NOT EXISTS ${OPTIONSFILE})

add_subdirectory("${BUCKETTOOLS_DIR}/ufc" "buckettools_ufc")

add_executable(${EXECUTABLE} "${PROJECT_SOURCE_DIR}/../../cpp/main.cpp")

set (EXTRA_LIBS ${EXTRA_LIBS} buckettools_cpp ${PREFIX}buckettools_ufc)

target_link_libraries(${EXECUTABLE} ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${EXTRA_LIBS} /usr/lib/libspud.a ${PYTHON_LIBRARIES})

add_custom_target(
  ${PREFIX}run
  COMMAND ${EXECUTABLE} ${OPTIONSFILE} 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${EXECUTABLE}
  )

set(RUNLIST ${RUNLIST} ${PREFIX}run)

add_custom_target(
  ${PREFIX}test
  COMMAND py.test ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )
add_dependencies(${PREFIX}test ${PREFIX}run)

