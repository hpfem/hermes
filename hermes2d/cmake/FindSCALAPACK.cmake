#
# ScaLAPACK and BLACS
#

FIND_LIBRARY(SCALAPACK_LIBRARY  NAMES scalapack scalapack-pvm scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-openmpi scalapack-lam 
                                PATHS /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

FIND_LIBRARY(BLACS_LIBRARY      NAMES blacs blacs-pvm blacs-mpi blacs-mpich blacs-mpich2 blacs-openmpi blacs-lam 
                                PATHS /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK DEFAULT_MSG SCALAPACK_LIBRARY BLACS_LIBRARY)

SET(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY} ${BLACS_LIBRARY})
