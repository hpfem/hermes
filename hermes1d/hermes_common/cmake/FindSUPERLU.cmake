#
# SUPERLU
#

FIND_PATH(SUPERLU_INCLUDE_DIR slu_ddefs.h ${SUPERLU_ROOT}/include /usr/include /usr/include/superlu /usr/local/include/superlu)

if(MSVC)
        set(SUPERLU_LIBRARY_NAME libsuperlu)
else(MSVC)
        set(SUPERLU_LIBRARY_NAME superlu)
endif(MSVC)
FIND_LIBRARY(SUPERLU_LIBRARY ${SUPERLU_LIBRARY_NAME} ${SUPERLU_ROOT}/lib /usr/lib /usr/local/lib/SUPERLU)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SUPERLU DEFAULT_MSG SUPERLU_INCLUDE_DIR SUPERLU_LIBRARY)
