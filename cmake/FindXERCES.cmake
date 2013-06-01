if(WIN32)
  set(XERCES_LIBRARY_NAMES xerces-c_3)
else(WIN32)
  set(XERCES_LIBRARY_NAMES xerces-c)
endif(WIN32)

FIND_PATH(XERCES_INCLUDE_DIR xercesc/sax/InputSource.hpp xercesc/dom/DOMDocument.hpp xercesc/dom/DOMErrorHandler.hpp ${XERCES_ROOT}/include)

if(64_BIT)
    FIND_LIBRARY(XERCES_LIBRARY ${XERCES_LIBRARY_NAMES} ${XERCES_ROOT}/lib/x64 /usr/lib64 /usr/local/lib64)
else(64_BIT) 
    FIND_LIBRARY(XERCES_LIBRARY ${XERCES_LIBRARY_NAMES} ${XERCES_ROOT}/lib /usr/lib /usr/local/lib)
endif(64_BIT)
    
# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XERCES DEFAULT_MSG XERCES_LIBRARY XERCES_INCLUDE_DIR)
