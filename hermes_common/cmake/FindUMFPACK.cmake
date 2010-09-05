#
# UMFPACK
#

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h ${UMFPACK_ROOT}/include /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse)
FIND_PATH(AMD_INCLUDE_DIR amd.h ${AMD_ROOT}/include /usr/include /usr/local/include/AMD /usr/include/suitesparse)

if(MSVC)
	set(UMFPACK_LIBRARY_NAME libumfpack)
	set(AMD_LIBRARY_NAME libamd)
else(MSVC)
	set(UMFPACK_LIBRARY_NAME umfpack)
	set(AMD_LIBRARY_NAME amd)
endif(MSVC)
FIND_LIBRARY(UMFPACK_LIBRARY ${UMFPACK_LIBRARY_NAME} ${UMFPACK_ROOT}/lib /usr/lib /usr/local/lib/UMFPACK)
FIND_LIBRARY(AMD_LIBRARY ${AMD_LIBRARY_NAME} ${AMD_ROOT}/lib /usr/lib /usr/local/lib/AMD)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK DEFAULT_MSG UMFPACK_INCLUDE_DIR AMD_INCLUDE_DIR UMFPACK_LIBRARY AMD_LIBRARY)
