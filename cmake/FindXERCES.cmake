#
# XERCES
#

FIND_PATH(XERCES_INCLUDE_DIR xercesc/sax/InputSource.hpp xercesc/dom/DOMDocument.hpp xercesc/dom/DOMErrorHandler.hpp ${XERCES_ROOT}/include)

if(WIN64)
    FIND_LIBRARY(XERCES_LIBRARY NAMES xerces-c_3 xerces-c PATHS ${XERCES_ROOT}/lib/x64 ${XERCES_ROOT}/lib)
else(WIN64) 
    FIND_LIBRARY(XERCES_LIBRARY NAMES xerces-c_3 xerces-c PATHS ${XERCES_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)
    
# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XERCES DEFAULT_MSG XERCES_LIBRARY XERCES_INCLUDE_DIR)
