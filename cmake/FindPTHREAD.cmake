#
# Pthread
#

if(MSVC)
	set(PTHREAD_LIBRARY_NAME pthreadVCE2)
else(MSVC)
	set(PTHREAD_LIBRARY_NAME pthread)
endif(MSVC)

FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${PTHREAD_ROOT}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PTHREAD DEFAULT_MSG PTHREAD_LIBRARY)
