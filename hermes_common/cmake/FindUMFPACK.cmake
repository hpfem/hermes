#
# UMFPACK
#

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h ${UMFPACK_ROOT}/include ${UMFPACK_ROOT}/Include /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse /opt/local/include/ufsparse)
FIND_PATH(AMD_INCLUDE_DIR amd.h ${AMD_ROOT}/include ${AMD_ROOT}/Include /usr/include /usr/local/include/AMD /usr/include/suitesparse /opt/local/include/ufsparse)

FIND_LIBRARY(UMFPACK_LIBRARY NAMES libumfpack umfpack PATHS ${UMFPACK_ROOT}/bin ${UMFPACK_ROOT}/Lib /usr/lib /usr/local/lib/UMFPACK) 
FIND_LIBRARY(AMD_LIBRARY NAMES libamd amd PATHS ${AMD_ROOT}/bin ${AMD_ROOT}/Lib /usr/lib /usr/local/lib/AMD)

SET(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
SET(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${AMD_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK DEFAULT_MSG UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS)

