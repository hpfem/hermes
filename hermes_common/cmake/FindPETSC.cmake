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

# You can specify your own version of the library instead of the one provided by
# Femhub by specifying the environment variables MY_PETSC_LIB_DIRS and 
# MY_PETSC_INC_DIRS.
IF ("$ENV{MY_PETSC_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_PETSC_INC_DIRS}" STREQUAL "")
  # When linking the library to stand-alone Hermes, you may also specify the 
  # variables directly in CMake.vars
  IF (NOT MY_PETSC_LIB_DIRS OR NOT MY_PETSC_INC_DIRS)
    # Alternatively, you may simply specify PETSC_ROOT in CMake.vars. This is 
    # the traditional way used also in the spkg files from the hpfem/solvers
    # repository and in the Hermes spkg.
    SET(MY_PETSC_LIB_DIRS ${PETSC_ROOT}/lib)
    SET(MY_PETSC_INC_DIRS ${PETSC_ROOT}/include)
  ENDIF (NOT MY_PETSC_LIB_DIRS OR NOT MY_PETSC_INC_DIRS)  
ELSE ("$ENV{MY_PETSC_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_PETSC_INC_DIRS}" STREQUAL "")
  SET(MY_PETSC_LIB_DIRS $ENV{MY_PETSC_LIB_DIRS})
  SET(MY_PETSC_INC_DIRS $ENV{MY_PETSC_INC_DIRS})
ENDIF ("$ENV{MY_PETSC_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_PETSC_INC_DIRS}" STREQUAL "")

# Try to find petsc.h in the root include directory.
FIND_PATH(COMMON_PETSC_INCLUDE_DIRS petsc.h PATHS ${MY_PETSC_INC_DIRS} NO_DEFAULT_PATH)
FIND_PATH(COMMON_PETSC_INCLUDE_DIRS petsc.h)

IF(NOT WITH_MPI AND EXISTS ${MY_PETSC_INC_DIRS}/mpiuni)
  IF(COMMON_PETSC_INCLUDE_DIRS)
    # Add path for the sequential emulation of MPI.
    SET(COMMON_PETSC_INCLUDE_DIRS ${COMMON_PETSC_INCLUDE_DIRS} ${MY_PETSC_INC_DIRS}/mpiuni)
  ELSE(COMMON_PETSC_INCLUDE_DIRS)
    # Set path for the sequential emulation of MPI.
    SET(COMMON_PETSC_INCLUDE_DIRS ${MY_PETSC_INC_DIRS}/mpiuni)
  ENDIF(COMMON_PETSC_INCLUDE_DIRS)
ENDIF(NOT WITH_MPI AND EXISTS ${MY_PETSC_INC_DIRS}/mpiuni)

SET(BASIC_PETSC_ARCH ${PETSC_ARCH})

IF(HERMES_COMMON_REAL) # Search for the real version of the library.

  # Look for the libraries only if they are not already in cache.
  IF(NOT PETSC_REAL_LIBRARIES)  
    SET(PETSC_ARCH ${BASIC_PETSC_ARCH}-real)

    # PETSc 3.1    
    SET(PETSC_DIR ${MY_PETSC_LIB_DIRS}/${PETSC_ARCH})
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
        MESSAGE( FATAL_ERROR
          "Real version of PETSC could not be found. Either disable it by setting 
           WITH_PETSC to NO in your CMake.vars file, or install it according to 
           instructions at\n
          <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/petsc.html>."
        )
      ENDIF (PETSC_FIND_REQUIRED)
    ENDIF (PETSC_FOUND) 
    
    # linux specific (?)
    SET(PETSC_REAL_LIBRARIES ${PETSC_REAL_LIBRARIES} dl 
      CACHE FILEPATH "PETSc libraries - real version")     
  ENDIF(NOT PETSC_REAL_LIBRARIES)
  
  # Export to cache.
  SET(PETSC_REAL_INCLUDE_DIRS ${PETSC_REAL_INCLUDE_DIRS} 
    CACHE PATH "PETSc include directories for the real version")
    
  UNSET(PETSC_LIB CACHE)  # PETSC_LIB don't needed any more - wipe out from cache.
ENDIF(HERMES_COMMON_REAL)
  
IF(HERMES_COMMON_COMPLEX)  # Search for the complex version of the library.  
  
  # Look for the libraries only if they are not already in cache.
  IF(NOT PETSC_CPLX_LIBRARIES)      
    # Reset the search flags.
    SET(PETSC_FOUND FALSE)
    
    SET(PETSC_ARCH ${BASIC_PETSC_ARCH}-complex)
    
    # PETSc 3.1    
    SET(PETSC_DIR ${MY_PETSC_LIB_DIRS}/${PETSC_ARCH})
    IF(EXISTS ${PETSC_DIR})
      FIND_LIBRARY(PETSC_LIB_C petsc ${PETSC_DIR}/lib NO_DEFAULT_PATH)
      FIND_LIBRARY(PETSC_LIB_C petsc)
                
      IF(COMMON_PETSC_INCLUDE_DIRS)
        # Add arch-specific include directory.
        SET(PETSC_CPLX_INCLUDE_DIRS ${PETSC_DIR}/include)
      ELSE(COMMON_PETSC_INCLUDE_DIRS)
        # petsc.h has not been found in the root include directory, search in the arch-specific one.
        FIND_PATH(PETSC_CPLX_INCLUDE_DIRS petsc.h PATHS ${PETSC_DIR}/include)
      ENDIF(COMMON_PETSC_INCLUDE_DIRS)
     	
      IF (PETSC_CPLX_INCLUDE_DIRS AND PETSC_LIB_C)
        SET(PETSC_FOUND TRUE)        
        # Set the complex version of the library (PETSc 3.1 is contained in a single libfile).
        SET(PETSC_CPLX_LIBRARIES ${PETSC_LIB_C})
      ENDIF (PETSC_CPLX_INCLUDE_DIRS AND PETSC_LIB_C)
      
    ENDIF(EXISTS ${PETSC_DIR})
    
    IF (PETSC_FOUND)    
      IF (NOT PETSC_FIND_QUIETLY)
		    MESSAGE(STATUS "Found complex version of PETSc: ${PETSC_DIR}")
      ENDIF (NOT PETSC_FIND_QUIETLY)
    ELSE (PETSC_FOUND)
      IF (PETSC_FIND_REQUIRED)
	      MESSAGE( FATAL_ERROR
          "Complex version of PETSC could not be found. Either disable it by setting 
           WITH_PETSC to NO in your CMake.vars file, or install it according to 
           instructions at\n
           <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/petsc.html>."
        )
      ENDIF (PETSC_FIND_REQUIRED)
    ENDIF (PETSC_FOUND) 
    
    # linux specific (?)
    SET(PETSC_CPLX_LIBRARIES ${PETSC_CPLX_LIBRARIES} dl CACHE FILEPATH "PETSc libraries - complex version")  
  ENDIF(NOT PETSC_CPLX_LIBRARIES)  

    # Export to cache.
  SET(PETSC_CPLX_INCLUDE_DIRS ${PETSC_CPLX_INCLUDE_DIRS} 
    CACHE PATH "PETSc include directories for the complex version")
    
  UNSET(PETSC_LIB_C CACHE)  # PETSC_LIB don't needed any more - wipe out from cache.
ENDIF(HERMES_COMMON_COMPLEX)
