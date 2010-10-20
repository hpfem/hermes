// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_SUPERLU
#include <superlu/slu_ddefs.h>

bool CommonSolverSuperLU::_solve(Matrix *mat, double *res)
{
    printf("SuperLU solver\n");

    int size = mat->get_size();
    int nnz = 0;

    int *Ap = NULL;
    int *Ai = NULL;
    double *Ax = NULL;

    CSCMatrix *Acsc = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acsc = new CSCMatrix(mcoo);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acsc = mcsc;
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acsc = new CSCMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    nnz = Acsc->get_nnz();
    Ap = Acsc->get_Ap();
    Ai = Acsc->get_Ai();
    Ax = Acsc->get_Ax();

    SuperMatrix A;
    SuperMatrix B;
    SuperMatrix L;      // factor L
    SuperMatrix U;      // factor U

    int nrhs, info;

    superlu_options_t options;
    SuperLUStat_t stat;

    // Set the default input options:
    /*
    options.Fact = DOFACT;
    options.Equil = YES;
    options.ColPerm = COLAMD;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
    options.SymmetricMode = NO;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;
    options.PrintStat = YES;
    */
    set_default_options(&options);

    // create csc matrix
    dCreate_CompCol_Matrix(&A, size, size, nnz, Ax, Ai, Ap,
                           SLU_NC, SLU_D, SLU_GE);
    // dPrint_CompCol_Matrix("A", &A);

    nrhs = 1;
    // create rhs matrix
    dCreate_Dense_Matrix(&B, size, nrhs, res, size,
                         SLU_DN, SLU_D, SLU_GE);
    // dPrint_Dense_Matrix("B", &B);

    // column permutation vector
    int *perm_c = intMalloc(size);
    // row permutations from partial pivoting
    int *perm_r = intMalloc(size);
    if (!perm_c) ABORT("Malloc fails for perm_c[].");
    if (!perm_r) ABORT("Malloc fails for perm_r[].");

    // initialize the statistics variables
    StatInit(&stat);
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    mem_usage_t mem_usage;
    if ( info == 0 )
    {
        // solution
        double *x = (double*) ((DNformat*) B.Store)->nzval;

        // copy result
        memcpy(res, x, size*sizeof(double));

        /*
        SCformat *Lstore = (SCformat *) L.Store;
        NCformat *Ustore = (NCformat *) U.Store;
        printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
        printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
        printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - mcsc.get_size());
        printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - mcsc.get_size())/nnz);

        mem_usage_t mem_usage;
        dQuerySpace(&L, &U, &mem_usage);
        printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        */
    }
    else
    {
        printf("dgssv() error returns INFO = %d\n", info);
        if (info <= size)
        {
            // factorization completes
            dQuerySpace(&L, &U, &mem_usage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        }
    }

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    // FIXME
    // Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    if (!dynamic_cast<CSCMatrix*>(mat))
        delete Acsc;
}

bool CommonSolverSuperLU::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSuperLU::solve(Matrix *mat, cplx *res) not implemented.");
}

#else

bool CommonSolverSuperLU::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverSuperLU::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSuperLU::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSuperLU::solve(Matrix *mat, cplx *res) not implemented.");
}

#endif
