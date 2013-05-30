#
# Pthread
#

if(MSVC)
	set(PTHREAD_LIBRARY_NAME pthreadVCE2)
else(MSVC)
	set(PTHREAD_LIBRARY_NAME pthread)
endif(MSVC)

if(64_BIT)
  FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${PTHREAD_ROOT}/lib/x64 /usr/lib64 /usr/local/lib64)
else(64_BIT)  
  FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${PTHREAD_ROOT}/lib /usr/lib /usr/local/lib)
endif(64_BIT)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PTHREAD DEFAULT_MSG PTHREAD_LIBRARY)
