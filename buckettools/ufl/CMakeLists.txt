
add_custom_target(
     generateufl
     python ${BucketTools_SOURCE_DIR}/scripts/ufl_from_options.py -s "temp" ${OPTIONSFILE}
     )
  
add_custom_command(
  OUTPUT SystemsWrapper.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_if_different ${BucketTools_BINARY_DIR}/ufl/SystemsWrapper.h.temp ${BucketTools_BINARY_DIR}/ufl/SystemsWrapper.h
  WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
  DEPENDS generateufl
  )

macro (compile_ufl ufl)
  add_custom_command(
    OUTPUT ${ufl}.ufl
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_if_different ${BucketTools_BINARY_DIR}/ufl/${ufl}.ufl.temp ${BucketTools_BINARY_DIR}/ufl/${ufl}.ufl
    WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
    DEPENDS generateufl
    )

  add_custom_command(
    OUTPUT ${ufl}.h ${ufl}.cpp
    COMMAND ffc 
    ARGS -l dolfin -O -f split ${ufl}.ufl
    WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
    DEPENDS generateufl
    )

  set (FFC_SOURCE ${FFC_SOURCE} ${ufl}.h)
  set (UFC_SOURCE ${UFC_SOURCE} ${ufl}.cpp)
endmacro (compile_ufl)

# We don't actually know what ufl files we've just generated so we need to find out what they're called!
execute_process (
  COMMAND ls *.ufl.temp
  WORKING_DIRECTORY ${BucketTools_BINARY_DIR}/ufl
  OUTPUT_VARIABLE UFL_SOURCE
  ERROR_VARIABLE UFL_ERR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_STRIP_TRAILING_WHITESPACE
  )

if(UFL_ERR)
  message(${UFL_ERR})
endif(UFL_ERR)

#foreach(f ${UFL_SOURCE})
#  compile_ufl (${f})
#endforeach(f)
compile_ufl (PoissonSimpleSolver)

add_custom_target(
  bucketffcufl
  SOURCES ${FFC_SOURCE}
  )

add_library(BucketUFCFiles ${UFC_SOURCE})
add_dependencies(BucketUFCFiles bucketffcufl generateufl SystemsWrapper.h)

