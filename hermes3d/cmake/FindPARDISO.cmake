#
# PARDISO
#
# set WITH_PARDISO to YES to enable pardiso support
# set PARDISO_LIB to the point to your pardiso library to use (use full path specification)
#
# If you are using parallel version of PARDISO, you need to say: set(WITH_OPENMP YES) in global
# CMake.vars. Depending on your configuration, you might need to link some additional libraries
# like gfortran (you can specify them in ADDITIONAL_LIB variable in your global CMake.vars file)
#

IF(EXISTS ${PARDISO_LIB})
	SET(PARDISO_LIBRARY ${PARDISO_LIB})
	SET(PARDISO_FOUND TRUE)
ENDIF(EXISTS ${PARDISO_LIB})

IF (PARDISO_FOUND)
	IF (NOT PARDISO_FIND_QUIETLY)
		MESSAGE(STATUS "Found PARDISO: ${PARDISO_LIBRARY}")
	ENDIF (NOT PARDISO_FIND_QUIETLY)
ELSE (PARDISO_FOUND)
	IF (PARDISO_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find PARDISO")
	ENDIF (PARDISO_FIND_REQUIRED)
ENDIF (PARDISO_FOUND)
