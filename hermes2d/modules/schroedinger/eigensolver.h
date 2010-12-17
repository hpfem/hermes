#ifndef __H2D_EIGENSOLVER_H
#define __H2D_EIGENSOLVER_H

#include "hermes2d.h"

class HERMES_API EigenSolver {
public:
    EigenSolver(const Matrix &A, const Matrix &B) {
        this->A = &A;
        this->B = &B;
    }

    void solve() {
    }

    void print_eigenvalues() {
    }

private:
    const Matrix *A;
    const Matrix *B;
};

#endif
