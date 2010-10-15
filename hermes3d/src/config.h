#define WITH_UMFPACK
/* #undef WITH_PARDISO */
/* #undef WITH_MUMPS */
/* #undef WITH_PETSC */
/* #undef WITH_HDF5 */
#define WITH_EXODUSII
/* #undef WITH_MPI */

// trilinos
#define WITH_TRILINOS
#define HAVE_AMESOS
#define HAVE_AZTECOO
#define HAVE_TEUCHOS
/* #undef HAVE_TEUCHOS_LINK */
/* #undef HAVE_TEUCHOS_BFD */
#define HAVE_EPETRA
#define HAVE_IFPACK
#define HAVE_ML
#define HAVE_NOX
#define HAVE_KOMPLEX

/* #undef TRACING */
#define DEBUG
/* #undef DEBUG_ORDER */

// elements
#define WITH_TETRA
#define WITH_HEX
/* #undef WITH_PRISM */

/* #undef PRELOADING */

/* --- */
#define PACKAGE_BUGREPORT "dandrs@unr.edu"
#define PACKAGE_NAME "Hermes3D"
#define PACKAGE_TARNAME "hermes3d"
#define PACKAGE_VERSION "0.0.1"
#define PACKAGE_STRING "Hermes3D v" PACKAGE_VERSION

