#
# PJLIB
#

FIND_PATH(PJLIB_INCLUDE_DIR pjlib.h ${PJLIB_ROOT}/include)

if(WIN64)
  FIND_LIBRARY(PJLIB_LIBRARY pjlib ${PJLIB_ROOT}/lib/x64 ${PJLIB_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(PJLIB_LIBRARY pjlib ${PJLIB_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PJLIB DEFAULT_MSG PJLIB_LIBRARY PJLIB_INCLUDE_DIR)
