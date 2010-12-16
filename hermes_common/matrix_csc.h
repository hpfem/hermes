#ifndef _CSC_MATRIX_H_
#define _CSC_MATRIX_H_

#include "solver/umfpack_solver.h"

class HERMES_API CSCMatrix: public UMFPackMatrix {
};

class HERMES_API AVector: public UMFPackVector {
};

#endif
