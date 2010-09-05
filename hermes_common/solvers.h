// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_SOLVERS_H
#define __HERMES_COMMON_SOLVERS_H

class Matrix;
class Vector;

// abstract class
class CommonSolver
{
public:
    virtual bool _solve(Matrix *mat, double *res) = 0;
    virtual bool _solve(Matrix *mat, cplx *res) = 0;
    virtual bool solve(Matrix *mat, Vector *res);
    inline char *get_log() { return log; }

private:
    char *log;
};

// c++ cg
class CommonSolverCG : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res)
    {
        _solve(mat, res, 1e-6, 1000);
    }
    bool _solve(Matrix *mat, double *res,
               double tol,
               int maxiter);
    bool _solve(Matrix *mat, cplx *res);
};
inline bool solve_linear_system_cg(Matrix *mat, double *res,
                                   double tolerance,
                                   int maxiter)
{
    CommonSolverCG solver;
    return solver._solve(mat, res, tolerance, maxiter);
}
inline bool solve_linear_system_cg(Matrix *mat, cplx *res)
{
    CommonSolverCG solver;
    return solver._solve(mat, res);
}

// c++ lu
class CommonSolverDenseLU : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_dense_lu(Matrix *mat, double *res)
{
    CommonSolverDenseLU solver;
    solver._solve(mat, res);
}
inline void solve_linear_system_dense_lu(Matrix *mat, cplx *res)
{
    CommonSolverDenseLU solver;
    solver._solve(mat, res);
}

// c++ umfpack - optional
class CommonSolverUmfpack : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_umfpack(Matrix *mat, double *res)
{
    CommonSolverUmfpack solver;
    solver._solve(mat, res);
}
inline void solve_linear_system_umfpack(Matrix *mat, cplx *res)
{
    CommonSolverUmfpack solver;
    solver._solve(mat, res);
}

// c++ sparselib - optional
class CommonSolverSparseLib : public CommonSolver
{
public:
    enum CommonSolverSparseLibSolver
    {
        HERMES_CommonSolverSparseLibSolver_ConjugateGradientSquared,
        CommonSolverSparseLibSolver_RichardsonIterativeRefinement
    };

    CommonSolverSparseLib()
    {
        tolerance = 1e-8;
        maxiter = 1000;
        method = HERMES_CommonSolverSparseLibSolver_ConjugateGradientSquared;
    }

    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
    inline void set_tolerance(double tolerance) { this->tolerance = tolerance; }
    inline void set_maxiter(int maxiter) { this->maxiter = maxiter; }
    inline void set_method(CommonSolverSparseLibSolver method) { this->method = method; }

private:
    double tolerance;
    int maxiter;
    CommonSolverSparseLibSolver method;
};
inline void solve_linear_system_sparselib_cgs(Matrix *mat, double *res, double tolerance = 1e-8, int maxiter = 1000)
{
    CommonSolverSparseLib solver;
    solver.set_method(CommonSolverSparseLib::HERMES_CommonSolverSparseLibSolver_ConjugateGradientSquared);
    solver.set_tolerance(tolerance);
    solver.set_maxiter(maxiter);
    solver._solve(mat, res);
}
inline void solve_linear_system_sparselib_ir(Matrix *mat, double *res, double tolerance = 1e-8, int maxiter = 1000)
{
    CommonSolverSparseLib solver;
    solver.set_method(CommonSolverSparseLib::CommonSolverSparseLibSolver_RichardsonIterativeRefinement);
    solver.set_tolerance(tolerance);
    solver.set_maxiter(maxiter);
    solver._solve(mat, res);
}

// c++ superlu - optional
class CommonSolverSuperLU : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve2(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_superlu(Matrix *mat, double *res)
{
    CommonSolverSuperLU solver;
    solver._solve(mat, res);
}

// python numpy - optional
class CommonSolverNumPy : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_numpy(Matrix *mat, double *res)
{
    CommonSolverNumPy solver;
    solver._solve(mat, res);
}
inline void solve_linear_system_numpy(Matrix *mat, cplx *res)
{
    CommonSolverNumPy solver;
    solver._solve(mat, res);
}

// python scipy - optional
class CommonSolverSciPyUmfpack : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_scipy_umfpack(Matrix *mat, double *res)
{
    CommonSolverSciPyUmfpack solver;
    solver._solve(mat, res);
}
inline void solve_linear_system_scipy_umfpack(Matrix *mat, cplx *res)
{
    CommonSolverSciPyUmfpack solver;
    solver._solve(mat, res);
}

class CommonSolverSciPyCG : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_scipy_cg(Matrix *mat, double *res)
{
    CommonSolverSciPyCG solver;
    solver._solve(mat, res);
}

class CommonSolverSciPyGMRES : public CommonSolver
{
public:
    bool _solve(Matrix *mat, double *res);
    bool _solve(Matrix *mat, cplx *res);
};
inline void solve_linear_system_scipy_gmres(Matrix *mat, double *res)
{
    CommonSolverSciPyGMRES solver;
    solver._solve(mat, res);
}

#endif
