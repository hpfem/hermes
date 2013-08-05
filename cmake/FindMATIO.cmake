#
# MATIO
# FROM http://sourceforge.net/projects/matio/
#

FIND_PATH(MATIO_INCLUDE_DIR matio.h ${MATIO_ROOT}/src ${MATIO_ROOT}/include /usr/local/include /usr/include)

if(WIN64)
  FIND_LIBRARY(MATIO_LIBRARY NAMES agros2d_3rd_party_matio matio libmatio PATHS ${MATIO_ROOT} ${MATIO_ROOT}/lib/x64)
else(WIN64)  
  FIND_LIBRARY(MATIO_LIBRARY NAMES agros2d_3rd_party_matio matio libmatio PATHS ${MATIO_ROOT} ${MATIO_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MATIO DEFAULT_MSG MATIO_LIBRARY MATIO_INCLUDE_DIR)
