#
# BSON
# FROM https://github.com/mongodb/mongo-c-driver
#

FIND_PATH(BSON_INCLUDE_DIR bson.h ${BSON_ROOT}/src ${BSON_ROOT}/include /usr/local/include /usr/include)

if(WIN64)
  FIND_LIBRARY(BSON_LIBRARY NAMES agros2d_3rd_party_bson bson PATHS ${BSON_ROOT} ${BSON_ROOT}/lib/x64)
else(WIN64)  
  FIND_LIBRARY(BSON_LIBRARY NAMES agros2d_3rd_party_bson bson PATHS ${BSON_ROOT} ${BSON_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BSON DEFAULT_MSG BSON_LIBRARY BSON_INCLUDE_DIR)
