if(WIN32)
  set(XERCES_LIBRARY_NAMES xerces-c_3)
else(WIN32)
  set(XERCES_LIBRARY_NAMES xerces-c)
endif(WIN32)

FIND_PATH(XERCES_INCLUDE_DIR xercesc/sax/InputSource.hpp xercesc/dom/DOMDocument.hpp xercesc/dom/DOMErrorHandler.hpp ${XERCES_ROOT}/include)

FIND_LIBRARY(XERCES_LIBRARY ${XERCES_LIBRARY_NAMES} ${XERCES_ROOT}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XERCES DEFAULT_MSG XERCES_LIBRARY XERCES_INCLUDE_DIR)
