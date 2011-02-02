#
# PETSc and MUMPS have different set of libraries for real and complex versions.
# PETSc also has specific include directories for real and complex versions, 
# apart from the architecture-independent ones. Following macros decide which 
# versions to use, according to which version of $HERMES we link the current 
# target $TRGT against. 
#
# Note that we cannot use 'include_directories' to add the architecture-specific
# include directories for PETSc as this function does not distinguish between
# build targets - these directories are thus excluded from dependency scanning, 
# but this is not an issue as long as there aren't any other libraries with
# same-named include files.
#
macro(SET_PETSC_FLAGS TRGT ARSP_INCLUDE_DIRS)
  if(WITH_PETSC)
    get_property(ARCH_SPECIFIC_FLAGS TARGET ${TRGT} PROPERTY COMPILE_FLAGS)
    foreach(_DIR ${${ARSP_INCLUDE_DIRS}})
      set(ARCH_SPECIFIC_FLAGS "${ARCH_SPECIFIC_FLAGS} -I${_DIR}")
    endforeach(_DIR)
    set_property(TARGET ${TRGT} PROPERTY COMPILE_FLAGS ${ARCH_SPECIFIC_FLAGS})
  endif(WITH_PETSC)
endmacro(SET_PETSC_FLAGS)

macro(PICK_REAL_OR_CPLX_LIBS HERMES TRGT)   
  if("${HERMES}" STREQUAL "${HERMES2D_REAL}" OR "${HERMES}" STREQUAL "${HERMES3D_REAL}")
    set(PETSC_LIBRARIES ${PETSC_REAL_LIBRARIES})
    set(MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
    SET_PETSC_FLAGS(${TRGT} PETSC_REAL_INCLUDE_DIRS)
  elseif("${HERMES}" STREQUAL "${HERMES2D_CPLX}" OR "${HERMES}" STREQUAL "${HERMES3D_CPLX}")
    set(PETSC_LIBRARIES ${PETSC_CPLX_LIBRARIES})
    set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})    
    SET_PETSC_FLAGS(${TRGT} PETSC_CPLX_INCLUDE_DIRS)
  endif("${HERMES}" STREQUAL "${HERMES2D_REAL}" OR "${HERMES}" STREQUAL "${HERMES3D_REAL}")
endmacro(PICK_REAL_OR_CPLX_LIBS)    
