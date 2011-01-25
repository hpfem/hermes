#include "eigensolver.h"

namespace Schroedinger {

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::ptr;
using Teuchos::null;
using Teuchos::rcp_dynamic_cast;

void wrap_CSC(const Ptr<Python> p, const std::string name,
        const Ptr<CSCMatrix> A)
{
    p->push_numpy_int_inplace("_IA", A->get_Ai(), A->get_nnz());
    p->push_numpy_int_inplace("_JA", A->get_Ap(), A->get_size()+1);
    p->push_numpy_double_inplace("_A", A->get_Ax(), A->get_nnz());
    p->push_int("n", A->get_size());
    p->exec("from scipy.sparse import csc_matrix\n");
    p->exec(name + " = csc_matrix((_A, _IA, _JA), shape=(n, n))");
}

PyMODINIT_FUNC initeigen(void); /*proto*/

void EigenSolver::solve() {
    // Support CSCMatrix only for now:
    RCP<CSCMatrix> A = rcp_dynamic_cast<CSCMatrix>(this->A);
    RCP<CSCMatrix> B = rcp_dynamic_cast<CSCMatrix>(this->B);
    Python p;
    wrap_CSC(ptr(&p), "A", A);
    wrap_CSC(ptr(&p), "B", B);
    initeigen();
    p.exec("from eigen import solve_eig_pysparse");

    printf("Solving the system A * x = lambda * B * x\n");
    p.exec("eigs = solve_eig_pysparse(A, B)");

    p.exec("energies = [E for E, eig in eigs]");
    printf("Energies:");
    p.exec("print energies");
}

} // namespace Schroedinger
