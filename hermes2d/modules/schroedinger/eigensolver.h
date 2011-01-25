#ifndef __H2D_EIGENSOLVER_H
#define __H2D_EIGENSOLVER_H

#include "hermes2d.h"

namespace Schroedinger {

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;

class HERMES_API EigenSolver {
public:
    EigenSolver(const RCP<Matrix> &A, const RCP<Matrix> &B);


    // Solves for 'n_eigs' eigenvectors, around the 'target_value'. Use
    // 'get_eigenvalue' and 'get_eigenvector' to retrieve the
    // eigenvalues/eigenvectors:
    void solve(int n_eigs=4, double target_value=-1);

    // Returns the number of calculated eigenvalues
    int get_n_eigs() {
        return this->n_eigs;
    }
    // Returns the i-th eigenvalue
    double get_eigenvalue(int i);
    // Returns the i-th eigenvector. A pointer will be returned into an
    // internal array, as well as the size of the vector. You don't own the
    // memory and it will be deallocated once the EigenSolver() class is
    // deleted. You need to make a copy of it if you want to store it
    // permanently.
    void get_eigenvector(int i, double **vec, int *n);

    void print_eigenvalues() {
        printf("Eigenvalues:\n");
        for (int i=0; i < this->get_n_eigs(); i++)
            printf("%3d: %f\n", i, this->get_eigenvalue(i));
    }

private:
    RCP<Matrix> A, B;
    int n_eigs;
    Python p;
};

}

#endif
