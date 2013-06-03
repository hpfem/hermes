#
# MUMPS
#
# set WITH_MUMPS to YES to enable MUMPS support
# set MUMPS_ROOT to point to the directory containing your MUMPS library
#

# You can specify your own version of the library instead of the one provided by
# Femhub by specifying the environment variables MY_MUMPS_LIB_DIRS and 
# MY_MUMPS_INC_DIRS.

# CMake maybe looks into the following paths by itself, but specifying them 
# explicitly doesn't hurt either.

IF(WIN32)
  MESSAGE(FATAL_ERROR "MUMPS only supported on Linux.")
ENDIF(WIN32)

SET(MUMPS_INCLUDE_SEARCH_PATH
	/usr/include
	/usr/include/mumps_seq
	/usr/local/include/
	/usr/local/include/mumps_seq
)

SET(MUMPS_LIB_SEARCH_PATH
	/usr/lib64
	/usr/lib
	/usr/local/lib/
)

FIND_PATH(MUMPS_INCLUDE_PATH  mumps_c_types.h ${MUMPS_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(MUMPS_MPISEQ_LIBRARY   mpiseq_seq        ${MUMPS_LIB_SEARCH_PATH})
FIND_LIBRARY(MUMPS_COMMON_LIBRARY   mumps_common_seq  ${MUMPS_LIB_SEARCH_PATH})
FIND_LIBRARY(MUMPS_PORD_LIBRARY     pord_seq          ${MUMPS_LIB_SEARCH_PATH})

FIND_PATH(MUMPS_MPISEQ_INCLUDE_PATH  mpi.h    ${MUMPS_INCLUDE_SEARCH_PATH})

SET(MUMPS_INCLUDE_PATH ${MUMPS_INCLUDE_PATH} ${MUMPS_MPISEQ_INCLUDE_PATH})

FIND_LIBRARY(MUMPSD_SEQ_LIBRARY dmumps_seq  ${MUMPS_LIB_SEARCH_PATH})
LIST(APPEND REQUIRED_REAL_LIBRARIES "MUMPSD_SEQ_LIBRARY")

FIND_LIBRARY(MUMPSZ_SEQ_LIBRARY zmumps_seq  ${MUMPS_LIB_SEARCH_PATH})
LIST(APPEND REQUIRED_CPLX_LIBRARIES "MUMPSZ_SEQ_LIBRARY")

LIST(APPEND REQUIRED_REAL_LIBRARIES "MUMPS_MPISEQ_LIBRARY")
LIST(APPEND REQUIRED_CPLX_LIBRARIES "MUMPS_MPISEQ_LIBRARY")  

LIST(APPEND REQUIRED_REAL_LIBRARIES "MUMPS_COMMON_LIBRARY" "MUMPS_PORD_LIBRARY")
LIST(APPEND REQUIRED_CPLX_LIBRARIES "MUMPS_COMMON_LIBRARY" "MUMPS_PORD_LIBRARY")

# Test if all the required libraries have been found. If they haven't, end with fatal error...
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(  MUMPS 
  "MUMPS could not be found. Either disable it by setting WITH_MUMPS to NO in
   your CMake.vars file, or install it according to instructions at\n
   <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/mumps.html>."
   ${REQUIRED_REAL_LIBRARIES} ${REQUIRED_CPLX_LIBRARIES} MUMPS_INCLUDE_PATH
) 

# ...if they have, append them all to the MUMPS_{REAL/CPLX}_LIBRARIES variable.
  FOREACH(_LIB ${REQUIRED_REAL_LIBRARIES})
    LIST(APPEND MUMPS_REAL_LIBRARIES ${${_LIB}})
  ENDFOREACH(_LIB ${REQUIRED_REAL_LIBRARIES})
  FOREACH(_LIB ${REQUIRED_CPLX_LIBRARIES})
    LIST(APPEND MUMPS_CPLX_LIBRARIES ${${_LIB}})
  ENDFOREACH(_LIB ${REQUIRED_CPLX_LIBRARIES})

# Finally, set MUMPS_INCLUDE_DIR to point to the MUMPS include directory.
SET(MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR} ${MUMPS_INCLUDE_PATH})
