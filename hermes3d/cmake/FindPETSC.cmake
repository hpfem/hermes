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
#   In PETSc source directory:
#     ./config/configure.py PETSC_ARCH=linux-cxx-real --with-clanguage=cxx
#     make PETSC_DIR=/opt/petsc/petsc-3.1-p4 PETSC_ARCH=linux-cxx-real all
#     ./config/configure.py PETSC_ARCH=linux-cxx-complex --with-clanguage=cxx --with-scalar-type=complex
#     make PETSC_DIR=/opt/petsc/petsc-3.1-p4 PETSC_ARCH=linux-mpicxx-complex all 
#
#   In hermes/hermes3d/CMake.vars:
#     set(WITH_PETSC YES)
#     set(PETSC_ROOT /opt/petsc/petsc-3.1-p4)
#     set(PETSC_ARCH linux-cxx)
#


# Try to find petsc.h in the root include directory.
FIND_PATH(COMMON_PETSC_INCLUDE_DIRS petsc.h PATHS ${PETSC_ROOT}/include)

IF(COMMON_PETSC_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_ROOT}/include/mpiuni)
  # Add path for the sequential emulation of MPI.
  SET(COMMON_PETSC_INCLUDE_DIRS ${COMMON_PETSC_INCLUDE_DIRS} ${PETSC_ROOT}/include/mpiuni)
ENDIF(COMMON_PETSC_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_ROOT}/include/mpiuni)

SET(BASIC_PETSC_ARCH ${PETSC_ARCH})

IF(H3D_REAL) # Search for the real version of the library.

  SET(PETSC_ARCH ${BASIC_PETSC_ARCH}-real)

  # PETSc 3.1    
  SET(PETSC_DIR ${PETSC_ROOT}/${PETSC_ARCH})
  IF(EXISTS ${PETSC_DIR})
    
    FIND_LIBRARY(PETSC_LIB petsc ${PETSC_DIR}/lib NO_DEFAULT_PATH)
    FIND_LIBRARY(PETSC_LIB petsc)
          
    IF(COMMON_PETSC_INCLUDE_DIRS)
      # Add arch-specific include directory.
      SET(PETSC_REAL_INCLUDE_DIRS ${PETSC_DIR}/include)
    ELSE(COMMON_PETSC_INCLUDE_DIRS)
      # petsc.h has not been found in the root include directory, search in the arch-specific one.
      FIND_PATH(PETSC_REAL_INCLUDE_DIRS petsc.h PATHS ${PETSC_DIR}/include)
    ENDIF(COMMON_PETSC_INCLUDE_DIRS)
   	  
    IF(PETSC_REAL_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_DIR}/include/mpiuni)
      # Add arch-specific path for mpiuni.
      SET(PETSC_REAL_INCLUDE_DIRS ${PETSC_REAL_INCLUDE_DIRS} ${PETSC_DIR}/include/mpiuni)
    ENDIF(PETSC_REAL_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_DIR}/include/mpiuni)
    
    IF (PETSC_REAL_INCLUDE_DIRS AND PETSC_LIB)
      SET(PETSC_FOUND TRUE)        
      # Set the real version of the library (PETSc 3.1 is contained in a single libfile).
      SET(PETSC_REAL_LIBRARIES ${PETSC_LIB})
    ENDIF (PETSC_REAL_INCLUDE_DIRS AND PETSC_LIB)
  
  ENDIF(EXISTS ${PETSC_DIR})  
      
  IF (PETSC_FOUND)    
    IF (NOT PETSC_FIND_QUIETLY)
		  MESSAGE(STATUS "Found real version of PETSc: ${PETSC_DIR}")
    ENDIF (NOT PETSC_FIND_QUIETLY)
  ELSE (PETSC_FOUND)
    IF (PETSC_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could not find real version of PETSc")
    ENDIF (PETSC_FIND_REQUIRED)
  ENDIF (PETSC_FOUND) 
  
  # linux specific (?)
  SET(PETSC_REAL_LIBRARIES ${PETSC_REAL_LIBRARIES} dl)     
ENDIF(H3D_REAL)
  
IF(H3D_COMPLEX)  # Search for the complex version of the library.  
  # Reset the search flags.
  SET(PETSC_FOUND FALSE)
  SET(PETSC_LIB PETSC_LIB-NOTFOUND)

  SET(PETSC_ARCH ${BASIC_PETSC_ARCH}-complex)
  
  # PETSc 3.1    
  SET(PETSC_DIR ${PETSC_ROOT}/${PETSC_ARCH})
  IF(EXISTS ${PETSC_DIR})
    
    FIND_LIBRARY(PETSC_LIB petsc ${PETSC_DIR}/lib NO_DEFAULT_PATH)
    FIND_LIBRARY(PETSC_LIB petsc)
          
    IF(COMMON_PETSC_INCLUDE_DIRS)
      # Add arch-specific include directory.
      SET(PETSC_CPLX_INCLUDE_DIRS ${PETSC_DIR}/include)
    ELSE(COMMON_PETSC_INCLUDE_DIRS)
      # petsc.h has not been found in the root include directory, search in the arch-specific one.
      FIND_PATH(PETSC_CPLX_INCLUDE_DIRS petsc.h PATHS ${PETSC_DIR}/include)
    ENDIF(COMMON_PETSC_INCLUDE_DIRS)
   	  
    IF(PETSC_CPLX_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_DIR}/include/mpiuni)
      # Add arch-specific path for mpiuni.
      SET(PETSC_CPLX_INCLUDE_DIRS ${PETSC_CPLX_INCLUDE_DIRS} ${PETSC_DIR}/include/mpiuni)
    ENDIF(PETSC_CPLX_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${PETSC_DIR}/include/mpiuni)
         		  
    IF (PETSC_CPLX_INCLUDE_DIRS AND PETSC_LIB)
      SET(PETSC_FOUND TRUE)        
      # Set the complex version of the library (PETSc 3.1 is contained in a single libfile).
      SET(PETSC_CPLX_LIBRARIES ${PETSC_LIB})
    ENDIF (PETSC_CPLX_INCLUDE_DIRS AND PETSC_LIB)
    
  ENDIF(EXISTS ${PETSC_DIR})
  
  IF (PETSC_FOUND)    
    IF (NOT PETSC_FIND_QUIETLY)
		  MESSAGE(STATUS "Found complex version of PETSc: ${PETSC_DIR}")
    ENDIF (NOT PETSC_FIND_QUIETLY)
  ELSE (PETSC_FOUND)
    IF (PETSC_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could not find complex version of PETSc")
    ENDIF (PETSC_FIND_REQUIRED)
  ENDIF (PETSC_FOUND) 
  
  # linux specific (?)
  SET(PETSC_CPLX_LIBRARIES ${PETSC_CPLX_LIBRARIES} dl)          
ENDIF(H3D_COMPLEX)



