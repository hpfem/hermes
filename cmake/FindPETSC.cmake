#
# PETSc
#
# When you configure and install PETSc, set PETSC_ROOT to some /root/dir/of/petsc/
# and PETSC_ARCH to petsc-arch-real (if you intend to solve real problems) and/or
# petsc-arch-complex (if you intend to solve complex problems). Then, in order to
# configure Hermes with PETSc, you need to say (in global CMake.vars):
#   set(WITH_PETSC YES)
#   set(PETSC_ROOT /root/dir/of/petsc/)
#   set(PETSC_ARCH petsc-arch)
#
# Example:
#   In PETSc source directory (this is automatically done by the
#   standalone-install script when using the prepackaged library from hpfem/solvers):
#     ./config/configure.py PETSC_ARCH=linux-cxx-real --with-clanguage=cxx
#     make PETSC_DIR=/opt/petsc/petsc-3.1-p7 PETSC_ARCH=linux-cxx-real all
#     ./config/configure.py PETSC_ARCH=linux-cxx-complex --with-clanguage=cxx --with-scalar-type=complex
#     make PETSC_DIR=/opt/petsc/petsc-3.1-p7 PETSC_ARCH=linux-mpicxx-complex all 
#
#   In hermes/CMake.vars:
#     set(WITH_PETSC YES)
#     set(PETSC_ROOT /opt/petsc/petsc-3.1-p7)
#     set(PETSC_ARCH linux-cxx)
#

IF(WIN32)
  MESSAGE(FATAL_ERROR "PETSc only supported on Linux.")
ENDIF(WIN32)

SET(COMMON_PETSC_INCLUDE_DIRS
    /usr/lib/petscdir/3.1/include
    /usr/include
    /usr/local/include
)

SET(COMMON_PETSC_LIB_DIRS
    /usr/lib/petscdir/3.1/linux-gnu-cxx-opt
    /usr/lib/petscdir/3.1/linux-gnu-c-opt/lib
    /usr/lib
    /usr/local/lib
    /usr/lib/petscdir/3.1/lib
)

FIND_PATH(PETSC_INCLUDE_DIRS petsc.h ${COMMON_PETSC_INCLUDE_DIRS})

# PETSc 3.1    
FIND_LIBRARY(PETSC_LIB_C petsc ${COMMON_PETSC_LIB_DIRS})

IF (PETSC_INCLUDE_DIRS AND PETSC_LIB_C)
  SET(PETSC_FOUND TRUE)        
  SET(PETSC_CPLX_LIBRARIES ${PETSC_LIB_C})
ENDIF (PETSC_INCLUDE_DIRS AND PETSC_LIB_C)


IF (PETSC_FOUND)    
IF (NOT PETSC_FIND_QUIETLY)
	    MESSAGE(STATUS "PETSc found: ${PETSC_LIB_C}")
ENDIF (NOT PETSC_FIND_QUIETLY)
ELSE (PETSC_FOUND)
IF (PETSC_FIND_REQUIRED)
      MESSAGE( FATAL_ERROR
  "PETSC could not be found. Either disable it by setting 
   WITH_PETSC to NO in your CMake.vars file, or install it according to 
   instructions at\n
   <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/petsc.html>."
)
ENDIF (PETSC_FIND_REQUIRED)
ENDIF (PETSC_FOUND)
