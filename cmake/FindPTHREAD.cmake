#
# Pthread
#

if(MSVC)
	set(PTHREAD_LIBRARY_NAME pthreadVCE2)
else(MSVC)
	set(PTHREAD_LIBRARY_NAME pthread)
endif(MSVC)

if(WIN64)
  FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${PTHREAD_ROOT}/lib/x64 ${PTHREAD_ROOT}/lib)
else(WIN64)  
  FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${PTHREAD_ROOT}/lib /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64)
endif(WIN64)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PTHREAD DEFAULT_MSG PTHREAD_LIBRARY)
