macro (compile_ufl ufl)
  add_custom_command(
    OUTPUT ${ufl}.ufl
    COMMAND python ${BucketTools_SOURCE_DIR}/scripts/complete_ufl_template.py
    ARGS ${OPTIONSFILE} ${BucketTools_SOURCE_DIR}/ufl/${ufl}.ufl.template ${BucketTools_BINARY_DIR}/ufl/${ufl}.ufl.temp
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_if_different ${BucketTools_BINARY_DIR}/ufl/${ufl}.ufl.temp ${BucketTools_BINARY_DIR}/ufl/${ufl}.ufl
    WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
    DEPENDS ${ufl}.ufl.template
    DEPENDS ${OPTIONSFILE}
    )

  add_custom_command(
    OUTPUT ${ufl}.h ${ufl}.cpp
    COMMAND ffc 
    ARGS -l dolfin -O -f split ${ufl}.ufl
    WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
    DEPENDS ${ufl}.ufl
    )

  set (FFC_SOURCE ${FFC_SOURCE} ${ufl}.h)
  set (UFC_SOURCE ${UFC_SOURCE} ${ufl}.cpp)
endmacro (compile_ufl)

compile_ufl (Divergence)

add_custom_target(
  bucketffcufl
  SOURCES ${FFC_SOURCE}
  )

add_library(BucketUFCFiles ${UFC_SOURCE})
add_dependencies(BucketUFCFiles bucketffcufl)

