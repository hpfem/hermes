#
# Cython
#

FIND_FILE(CYTHON_BIN cython /usr/bin /usr/local/bin) 

IF (CYTHON_BIN)
	SET(Cython_FOUND TRUE)
ENDIF (CYTHON_BIN)


IF (Cython_FOUND)
	IF (NOT Cython_FIND_QUIETLY)
		MESSAGE(STATUS "Found cython: ${CYTHON_BIN}")
	ENDIF (NOT Cython_FIND_QUIETLY)
ELSE (Cython_FOUND)
	IF (Cython_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find cython")
	ENDIF (Cython_FIND_REQUIRED)
ENDIF (Cython_FOUND)

