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

    printf("Solving the system A * x = lambda * B * x\n");
}

} // namespace Schroedinger
