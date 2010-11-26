#
# SuperLU
#

IF(SUPERLU_MT AND WITH_OPENMP)
  SET(POST _mt_OPENMP)
ELSEIF(SUPERLU_MT)
  SET(POST _mt_PTHREAD)
ENDIF(SUPERLU_MT AND WITH_OPENMP)

IF(POST)
  FIND_PATH(SUPERLU_INCLUDE_DIR pdsp_defs.h ${SUPERLU_ROOT}/include/superlu_mt NO_DEFAULT_PATH)
  FIND_PATH(SUPERLU_INCLUDE_DIR pdsp_defs.h /usr/include/superlu_mt /usr/local/include/superlu_mt)
ELSE(POST)
  FIND_PATH(SUPERLU_INCLUDE_DIR slu_ddefs.h ${SUPERLU_ROOT}/include/superlu NO_DEFAULT_PATH)
  FIND_PATH(SUPERLU_INCLUDE_DIR slu_ddefs.h /usr/include /usr/include/superlu /usr/local/include/superlu)
ENDIF(POST)

IF(MSVC)
  SET(PRE lib)
ENDIF(MSVC)

FIND_LIBRARY( SUPERLU_LIBRARY ${PRE}superlu${POST} ${SUPERLU_ROOT}/lib NO_DEFAULT_PATH)
FIND_LIBRARY( SUPERLU_LIBRARY ${PRE}superlu${POST} /usr/lib /usr/lib/superlu /usr/local/lib/superlu)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SUPERLU DEFAULT_MSG SUPERLU_LIBRARY SUPERLU_INCLUDE_DIR)
