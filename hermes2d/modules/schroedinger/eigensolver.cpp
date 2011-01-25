#include "eigensolver.h"

namespace Schroedinger {

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;
using Teuchos::rcp_dynamic_cast;

void EigenSolver::solve() {
    // Support CSCMatrix only for now:
    RCP<CSCMatrix> A = rcp_dynamic_cast<CSCMatrix>(this->A);
    RCP<CSCMatrix> B = rcp_dynamic_cast<CSCMatrix>(this->B);
    Python p;
    p.push_numpy_int_inplace("Ap", A->get_Ap(), A->get_size()+1);
    p.push_numpy_int_inplace("Ai", A->get_Ai(), A->get_nnz());
    p.push_numpy_double_inplace("Ax", A->get_Ax(), A->get_nnz());

    printf("Solving the system A * x = lambda * B * x\n");
}

} // namespace Schroedinger
