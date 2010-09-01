#ifndef __H2D_DUMMY_SOLVER_H
#define __H2D_DUMMY_SOLVER_H

#include "hermes2d.h"
#include "solver.h"
#include "config.h"

#if USE_UMFPACK
#include "solver_umfpack.h"
class DummySolver: public UmfpackSolver
{
};
#else
class DummySolver: public Solver
{
    virtual bool is_row_oriented() {
        return false;
    }
    virtual bool handles_symmetry() {
        return false;
    }
    virtual bool solve(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool
            sym, scalar* RHS, scalar* vec) {
        return true;
    }
};
#endif

#endif
