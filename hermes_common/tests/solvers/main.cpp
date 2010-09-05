#include <iostream>
#include <stdexcept>

#include "matrix.h"
#include "solvers.h"

#define EPS 1e-12

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_solver_dense_lu1()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_dense_lu(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver_dense_lu2()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_dense_lu(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    DenseMatrix B(&A);
    for (int i=0; i < 4; i++) res[i] = 1.;

    solve_linear_system_dense_lu(&B, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver_cg()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    _assert(solve_linear_system_cg(&A, res, EPS, 2));
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver_scipy_1()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    for (int i=0; i < 4; i++) res[i] = 1.;
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    for (int i=0; i < 4; i++) res[i] = 1.;
    solve_linear_system_scipy_cg(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver_scipy_2()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);

    res[0] = 8.;
    res[1] = 45.;
    res[2] = -3.;
    res[3] = 3.;
    res[4] = 19.;
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);

    res[0] = 8.;
    res[1] = 45.;
    res[2] = -3.;
    res[3] = 3.;
    res[4] = 19.;
    solve_linear_system_scipy_gmres(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);
}

void test_solver_scipy_3()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    cplx res[5];

    res[0] = cplx(8.);
    res[1] = cplx(45.);
    res[2] = cplx(-3.);
    res[3] = cplx(3.);
    res[4] = cplx(19.);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 1.) < EPS);
    _assert(fabs(res[1].real() - 2.) < EPS);
    _assert(fabs(res[2].real() - 3.) < EPS);
    _assert(fabs(res[3].real() - 4.) < EPS);
    _assert(fabs(res[4].real() - 5.) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);
    _assert(fabs(res[2].imag() - 0.) < EPS);
    _assert(fabs(res[3].imag() - 0.) < EPS);
    _assert(fabs(res[4].imag() - 0.) < EPS);

    res[0] = cplx(8.);
    res[1] = cplx(45.);
    res[2] = cplx(-3.);
    res[3] = cplx(3.);
    res[4] = cplx(19.);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 1.) < EPS);
    _assert(fabs(res[1].real() - 2.) < EPS);
    _assert(fabs(res[2].real() - 3.) < EPS);
    _assert(fabs(res[3].real() - 4.) < EPS);
    _assert(fabs(res[4].real() - 5.) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);
    _assert(fabs(res[2].imag() - 0.) < EPS);
    _assert(fabs(res[3].imag() - 0.) < EPS);
    _assert(fabs(res[4].imag() - 0.) < EPS);
}

void test_solver_scipy_4()
{
    CooMatrix A(2, true);
    A.add(0, 0, cplx(1, 1));
    A.add(0, 1, cplx(2, 2));
    A.add(1, 0, cplx(3, 3));
    A.add(1, 1, cplx(4, 4));

    cplx res[2];

    //----------------------

    res[0] = cplx(1);
    res[1] = cplx(2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.25) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - (-0.25)) < EPS);

    res[0] = cplx(1);
    res[1] = cplx(2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.25) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - (-0.25)) < EPS);

    //----------------------

    res[0] = cplx(1, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.5) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);

    res[0] = cplx(1, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.5) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);

    //----------------------

    res[0] = cplx(2, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - (-1)) < EPS);
    _assert(fabs(res[1].real() - 1.25) < EPS);
    _assert(fabs(res[0].imag() - 1.) < EPS);
    _assert(fabs(res[1].imag() - (-0.75)) < EPS);

    res[0] = cplx(2, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - (-1)) < EPS);
    _assert(fabs(res[1].real() - 1.25) < EPS);
    _assert(fabs(res[0].imag() - 1.) < EPS);
    _assert(fabs(res[1].imag() - (-0.75)) < EPS);
}

void test_solver_umfpack_real()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_umfpack(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);
}

void test_solver_umfpack_imag()
{
    CooMatrix A(2, true);
    A.add(0, 0, cplx(1, 1));
    A.add(0, 1, cplx(2, 2));
    A.add(1, 0, cplx(3, 3));
    A.add(1, 1, cplx(4, 4));

    cplx res[2];

    //----------------------

    res[0] = cplx(1);
    res[1] = cplx(2);
    solve_linear_system_umfpack(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.25) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - (-0.25)) < EPS);
}

void test_solver_sparselib_cgs()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_sparselib_cgs(&A, res, 1e-14);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);
}

void test_solver_sparselib_ir()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 1, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_sparselib_ir(&A, res, 1e-14);
    _assert(fabs(res[0] - 1.24489795918367) < EPS);
    _assert(fabs(res[1] - 1.83673469387755) < EPS);
    _assert(fabs(res[2] - 3.00000000000000) < EPS);
    _assert(fabs(res[3] - 3.91836734693878) < EPS);
    _assert(fabs(res[4] - 5.65306122448980) < EPS);
}

void test_solver_superlu()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_superlu(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);
}

int main(int argc, char* argv[])
{
    try {
        // SparseLib++
        test_solver_sparselib_cgs();
        test_solver_sparselib_ir();

        // Hermes Common
        test_solver_dense_lu1();
        test_solver_dense_lu2();
        test_solver_cg();

        // NumPy + SciPy
#ifdef COMMON_WITH_SCIPY
        test_solver_scipy_1();
        test_solver_scipy_2();
        test_solver_scipy_3();
        test_solver_scipy_4();
#endif
        // UMFPACK
#ifdef COMMON_WITH_UMFPACK
        test_solver_umfpack_real();
        test_solver_umfpack_imag();
#endif

        // SuperLU
#ifdef COMMON_WITH_SUPERLU
        test_solver_superlu();
#endif

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
