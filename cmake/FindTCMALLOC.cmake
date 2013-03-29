#
# TCMALLOC
# FROM http://code.google.com/p/google-perftools/
#

SET(MY_TCMALLOC_LIB_DIRS ${TCMALLOC_ROOT}/lib) 
SET(MY_TCMALLOC_INC_DIRS ${TCMALLOC_ROOT}/include)

FIND_PATH(TCMALLOC_INCLUDE_DIR tcmalloc.h ${MY_TCMALLOC_INC_DIRS} /usr/local/include/google /usr/include/google)

FIND_LIBRARY(TCMALLOC_LIBRARY  NAMES tcmalloc_minimal libtcmalloc_minimal PATHS ${MY_TCMALLOC_LIB_DIRS} /usr/lib /usr/local/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  TCMALLOC "TCMALLOC could not be found." 
   TCMALLOC_LIBRARY TCMALLOC_INCLUDE_DIR
)