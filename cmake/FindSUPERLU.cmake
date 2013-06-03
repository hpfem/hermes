#
# SuperLU
#

# You can specify your own version of the library instead of the one provided by
# Femhub by specifying the environment variables MY_SUPERLU_LIB_DIRS and 
# MY_SUPERLU_INC_DIRS.
IF ("$ENV{MY_SUPERLU_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_SUPERLU_INC_DIRS}" STREQUAL "")
  # When linking the library to stand-alone Hermes, you may also specify the 
  # variables directly in CMake.vars
  IF (NOT MY_SUPERLU_LIB_DIRS OR NOT MY_SUPERLU_INC_DIRS)
    # Alternatively, you may simply specify SUPERLU_ROOT in CMake.vars. This is 
    # the traditional way used also in the spkg files from the hpfem/solvers
    # repository and in the Hermes spkg.
    IF(WIN64)
      SET(MY_SUPERLU_LIB_DIRS ${SUPERLU_ROOT}/lib/x64 ${SUPERLU_ROOT}/lib)
    ELSE(WIN64)
      SET(MY_SUPERLU_LIB_DIRS ${SUPERLU_ROOT}/lib)
    ENDIF(WIN64)
    SET(MY_SUPERLU_INC_DIRS ${SUPERLU_ROOT}/include)
  ENDIF (NOT MY_SUPERLU_LIB_DIRS OR NOT MY_SUPERLU_INC_DIRS)
ELSE ("$ENV{MY_SUPERLU_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_SUPERLU_INC_DIRS}" STREQUAL "")
  SET(MY_SUPERLU_LIB_DIRS $ENV{MY_SUPERLU_LIB_DIRS})
  SET(MY_SUPERLU_INC_DIRS $ENV{MY_SUPERLU_INC_DIRS})
ENDIF ("$ENV{MY_SUPERLU_LIB_DIRS}" STREQUAL "" OR "$ENV{MY_SUPERLU_INC_DIRS}" STREQUAL "") 

IF(SUPERLU_MT AND WITH_OPENMP)
  SET(POST _mt_OPENMP)
ELSEIF(SUPERLU_MT)
  SET(POST _mt_PTHREAD)
ENDIF(SUPERLU_MT AND WITH_OPENMP)

IF(POST)
  FIND_PATH(SUPERLU_INCLUDE_DIR pdsp_defs.h ${MY_SUPERLU_INC_DIRS} NO_DEFAULT_PATH)
  FIND_PATH(SUPERLU_INCLUDE_DIR pdsp_defs.h /usr/include/superlu_mt /usr/local/include/superlu_mt)
ELSE(POST)
  FIND_PATH(SUPERLU_INCLUDE_DIR slu_ddefs.h ${MY_SUPERLU_INC_DIRS} NO_DEFAULT_PATH)
  FIND_PATH(SUPERLU_INCLUDE_DIR slu_ddefs.h /usr/include /usr/include/superlu /usr/local/include/superlu)
ENDIF(POST)

IF(MSVC)
  SET(PRE lib)
ENDIF(MSVC)

FIND_LIBRARY( SUPERLU_LIBRARY ${PRE}superlu${POST} ${MY_SUPERLU_LIB_DIRS} NO_DEFAULT_PATH)
FIND_LIBRARY( SUPERLU_LIBRARY ${PRE}superlu${POST} /usr/lib /usr/lib/superlu /usr/local/lib/superlu)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  SUPERLU 
   "SUPERLU could not be found. Please install it according to instructions at\n
   < http://hpfem.org/hermes/doc/src/installation/matrix_solvers/superlu.html >\n
   and/or provide path to its root directory by setting variable SUPERLU_ROOT 
   in the CMake.vars file."
  SUPERLU_LIBRARY SUPERLU_INCLUDE_DIR
)
