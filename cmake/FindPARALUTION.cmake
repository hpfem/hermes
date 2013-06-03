#
# PARALUTION
# FROM http://www.paralution.com/
#

FIND_PATH(PARALUTION_INCLUDE_DIR paralution.hpp ${PARALUTION_ROOT}/include ${PARALUTION_ROOT}/inc /usr/local/include/google /usr/include/google)

if(WIN64)
  FIND_LIBRARY(PARALUTION_LIBRARY paralution ${PARALUTION_ROOT}/lib/x64 ${PARALUTION_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(PARALUTION_LIBRARY paralution  ${PARALUTION_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PARALUTION DEFAULT_MSG PARALUTION_LIBRARY PARALUTION_INCLUDE_DIR)
