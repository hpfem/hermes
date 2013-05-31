#ifndef __HERMES_CONFIG_H_
#define __HERMES_CONFIG_H_

// OpenMP
#cmakedefine NUM_THREADS ${NUM_THREADS}
#cmakedefine WITH_OPENMP

#cmakedefine HAVE_FMEMOPEN
#cmakedefine HAVE_LOG2
#cmakedefine EXTREME_QUAD

#cmakedefine WITH_UMFPACK
#cmakedefine WITH_PARALUTION
#cmakedefine WITH_MUMPS
#cmakedefine WITH_SUPERLU
#cmakedefine WITH_PETSC
#cmakedefine WITH_HDF5
#cmakedefine WITH_EXODUSII
#cmakedefine WITH_MPI

#cmakedefine WITH_TC_MALLOC
#cmakedefine UMFPACK_LONG_INT

// stacktrace
#cmakedefine WITH_STACKTRACE
#cmakedefine EXECINFO_FOUND

// trilinos
#cmakedefine WITH_TRILINOS
#cmakedefine HAVE_AMESOS
#cmakedefine HAVE_AZTECOO
#cmakedefine HAVE_TEUCHOS
#cmakedefine HAVE_EPETRA
#cmakedefine HAVE_IFPACK
#cmakedefine HAVE_ML
#cmakedefine HAVE_NOX
#cmakedefine HAVE_KOMPLEX

// no logo
#cmakedefine HERMES_NO_LOGO

// GLUT
#cmakedefine NOGLUT

#endif