#ifndef __H2D_EIGENSOLVER_H
#define __H2D_EIGENSOLVER_H

#include "hermes2d.h"

class HERMES_API EigenSolver {
public:
    EigenSolver(const RCP<Matrix> &A, const RCP<Matrix> &B) {
        this->A = A;
        this->B = B;
    }

    void solve() {
        printf("Solving the system A * x = lambda * B * x\n");
    }

    void print_eigenvalues() {
    }

private:
    RCP<Matrix> A, B;
};

#endif
