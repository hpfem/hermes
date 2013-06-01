#
# PARALUTION
# FROM http://www.paralution.com/
#

FIND_PATH(PARALUTION_INCLUDE_DIR paralution.hpp ${PARALUTION_ROOT}/include ${PARALUTION_ROOT}/inc /usr/local/include/google /usr/include/google)
message(${PARALUTION_ROOT})
if(64_BIT)
  FIND_LIBRARY(PARALUTION_LIBRARY paralution ${PARALUTION_ROOT}/lib/x64 /usr/lib64 /usr/local/lib64)
else(64_BIT)  
  FIND_LIBRARY(PARALUTION_LIBRARY paralution  ${PARALUTION_ROOT}/lib /usr/lib /usr/local/lib)
endif(64_BIT)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PARALUTION DEFAULT_MSG PARALUTION_LIBRARY TCMALLOC_INCLUDE_DIR)
