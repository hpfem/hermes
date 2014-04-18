#
# LIBIBERTY
#

FIND_PATH(LIBIBERTY_INCLUDE_DIR libiberty.h ${LIBIBERTY_ROOT}/include/libiberty /usr/include/libiberty /usr/local/include/libiberty)

if(WIN64)
    FIND_LIBRARY(LIBIBERTY_LIBRARY NAMES libiberty PATHS ${LIBIBERTY_ROOT}/lib/x64 ${LIBIBERTY_ROOT}/lib)
else(WIN64) 
    FIND_LIBRARY(LIBIBERTY_LIBRARY NAMES libiberty.a PATHS ${LIBIBERTY_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64 /usr/lib/x86_64-linux-gnu /usr/lib/x86-linux-gnu)
endif(WIN64)
    
# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBIBERTY DEFAULT_MSG LIBIBERTY_LIBRARY LIBIBERTY_INCLUDE_DIR)
