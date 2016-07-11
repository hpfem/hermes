#
# XERCES
#

FIND_PATH(XERCES_INCLUDE_DIR xercesc/sax/InputSource.hpp xercesc/dom/DOMDocument.hpp xercesc/dom/DOMErrorHandler.hpp ${XERCES_ROOT}/include)

if(WIN64)
    FIND_LIBRARY(XERCES_LIBRARY NAMES xerces-c_3 xerces-c PATHS ${XERCES_ROOT}/lib/x64 ${XERCES_ROOT}/lib/x64/vc-12.0 ${XERCES_ROOT}/lib/x64/vc-11.0 ${XERCES_ROOT}/lib/x64/vc-10.0 ${XERCES_ROOT}/lib ${XERCES_ROOT}/lib/vc-12.0 ${XERCES_ROOT}/lib/vc-11.0 ${XERCES_ROOT}/lib/vc-10.0)
else(WIN64) 
    FIND_LIBRARY(XERCES_LIBRARY NAMES xerces-c_3 xerces-c PATHS ${XERCES_ROOT}/lib ${XERCES_ROOT}/lib/vc-12.0 ${XERCES_ROOT}/lib/vc-11.0 ${XERCES_ROOT}/lib/vc-10.0 /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)
    
# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XERCES DEFAULT_MSG XERCES_LIBRARY XERCES_INCLUDE_DIR)
