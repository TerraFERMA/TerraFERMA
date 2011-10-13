set(EXECUTABLE "diamond_validation")

set(PREFIX ${EXECUTABLE})

project($EXECUTABLE)

add_custom_target(
  ${PREFIX}_run
  COMMAND python ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )

add_custom_target(
  ${PREFIX}_test
  COMMAND py.test -vs ${PROJECT_SOURCE_DIR}/tests.py 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests.py
  )

