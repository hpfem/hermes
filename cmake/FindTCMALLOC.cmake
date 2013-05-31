#
# TCMALLOC
# FROM http://code.google.com/p/google-perftools/
#

FIND_PATH(TCMALLOC_INCLUDE_DIR tcmalloc.h ${TCMALLOC_ROOT}/include /usr/local/include/google /usr/include/google)

if(64_BIT)
  FIND_LIBRARY(TCMALLOC_LIBRARY  NAMES tcmalloc_minimal libtcmalloc_minimal PATHS ${TCMALLOC_ROOT}/lib/x64 /usr/lib64 /usr/local/lib64)
else(64_BIT)  
  FIND_LIBRARY(TCMALLOC_LIBRARY  NAMES tcmalloc_minimal libtcmalloc_minimal PATHS ${TCMALLOC_ROOT}/lib /usr/lib /usr/local/lib)
endif(64_BIT)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  TCMALLOC DEFAULT_MSG TCMALLOC_LIBRARY TCMALLOC_INCLUDE_DIR)