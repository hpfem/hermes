#
# TCMALLOC
# FROM http://code.google.com/p/google-perftools/
#

FIND_PATH(TCMALLOC_INCLUDE_DIR tcmalloc.h ${TCMALLOC_ROOT}/include /usr/local/include/google /usr/include/google)

if(WIN64)
  FIND_LIBRARY(TCMALLOC_LIBRARY  NAMES tcmalloc_minimal libtcmalloc_minimal PATHS ${TCMALLOC_ROOT}/lib/x64 ${TCMALLOC_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(TCMALLOC_LIBRARY  NAMES tcmalloc_minimal libtcmalloc_minimal PATHS ${TCMALLOC_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  TCMALLOC DEFAULT_MSG TCMALLOC_LIBRARY TCMALLOC_INCLUDE_DIR)