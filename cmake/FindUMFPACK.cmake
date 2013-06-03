#
# UMFPACK
#

# You can specify your own version of the library instead of the one provided by
# Femhub by specifying the environment variables MY_UMFPACK_LIB_DIRS and 
# MY_UMFPACK_INC_DIRS.
IF ("$ENV{MY_UMFPACK_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_UMFPACK_INC_DIRS}" STREQUAL "")
  # When linking the library to stand-alone Hermes, you may also specify the 
  # variables directly in CMake.vars
  IF (NOT MY_UMFPACK_LIB_DIRS OR NOT MY_UMFPACK_INC_DIRS)
    # Alternatively, you may simply specify UMFPACK_ROOT in CMake.vars. This is 
    # the traditional way used also in the spkg files from the hpfem/solvers
    # repository and in the Hermes spkg.
    if(WIN64)
        SET(MY_UMFPACK_LIB_DIRS ${UMFPACK_ROOT}/lib/x64 ${UMFPACK_ROOT}/lib) 
    else(WIN64) 
        SET(MY_UMFPACK_LIB_DIRS ${UMFPACK_ROOT}/lib) 
    endif(WIN64)
    SET(MY_UMFPACK_INC_DIRS ${UMFPACK_ROOT}/include)
  ENDIF (NOT MY_UMFPACK_LIB_DIRS OR NOT MY_UMFPACK_INC_DIRS)
ELSE ("$ENV{MY_UMFPACK_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_UMFPACK_INC_DIRS}" STREQUAL "")
  SET(MY_UMFPACK_LIB_DIRS $ENV{MY_UMFPACK_LIB_DIRS})
  SET(MY_UMFPACK_INC_DIRS $ENV{MY_UMFPACK_INC_DIRS})
ENDIF ("$ENV{MY_UMFPACK_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_UMFPACK_INC_DIRS}" STREQUAL "")

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h ${MY_UMFPACK_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(AMD_INCLUDE_DIR     amd.h     ${MY_UMFPACK_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse /opt/local/include/ufsparse)
FIND_PATH(AMD_INCLUDE_DIR     amd.h     /usr/include /usr/local/include/AMD /usr/include/suitesparse /opt/local/include/ufsparse)

FIND_LIBRARY(UMFPACK_LIBRARY  NAMES libumfpack umfpack  PATHS ${MY_UMFPACK_LIB_DIRS}  NO_DEFAULT_PATH) 
FIND_LIBRARY(SSC_LIBRARY  NAMES libsuitesparseconfig suitesparseconfig  PATHS ${MY_UMFPACK_LIB_DIRS}  NO_DEFAULT_PATH) 
FIND_LIBRARY(AMD_LIBRARY      NAMES libamd amd          PATHS ${MY_UMFPACK_LIB_DIRS}  NO_DEFAULT_PATH)

FIND_LIBRARY(UMFPACK_LIBRARY  NAMES libumfpack umfpack  PATHS /usr/lib /usr/local/lib/UMFPACK)
FIND_LIBRARY(SSC_LIBRARY  NAMES libsuitesparseconfig suitesparseconfig  PATHS /usr/lib /usr/local/lib/UMFPACK)
FIND_LIBRARY(AMD_LIBRARY      NAMES libamd amd          PATHS /usr/lib /usr/local/lib/AMD)

if(${SSC_LIBRARY} STREQUAL "SSC_LIBRARY-NOTFOUND")
  set(SSC_LIBRARY "")
endif(${SSC_LIBRARY} STREQUAL "SSC_LIBRARY-NOTFOUND")

SET(UMFPACK_INCLUDE_DIRS  ${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
SET(UMFPACK_LIBRARIES     ${UMFPACK_LIBRARY} ${SSC_LIBRARY} ${AMD_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  UMFPACK
   "UMFPACK could not be found. Please install it according to instructions at\n
   < http://hpfem.org/hermes/doc/src/installation/matrix_solvers/umfpack.html >\n
   and/or provide path to its root directory by setting variable UMFPACK_ROOT 
   in the CMake.vars file." 
   UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS
)