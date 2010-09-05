#include <iostream>
#include <stdexcept>

#include "python_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

Python *p;

void test_basic()
{
    p->push("i", c2py_int(5));
    p->exec("i = i*2");
    int i = py2c_int(p->pull("i"));
    _assert(i == 10);
}

void test_numpy()
{
    // double arrays
    p->exec("from numpy import array");
    p->exec("a = array(range(20), dtype='double')");
    p->exec("assert a.strides == (8,)");
    p->exec("b = a[::5]");
    p->exec("assert b.strides == (40,)");
    double *A;
    int n;
    numpy2c_double_inplace(p->pull("a"), &A, &n);
    _assert(n == 20);
    _assert(
            (fabs(A[0] - 0.)  < 1e-10) &&
            (fabs(A[1] - 1.)  < 1e-10) &&
            (fabs(A[2] - 2.) < 1e-10) &&
            (fabs(A[3] - 3.) < 1e-10)
           );
    numpy2c_double_inplace(p->pull("b"), &A, &n);
    _assert(n == 4);
    _assert(
            (fabs(A[0] - 0.)  < 1e-10) &&
            (fabs(A[1] - 5.)  < 1e-10) &&
            (fabs(A[2] - 10.) < 1e-10) &&
            (fabs(A[3] - 15.) < 1e-10)
           );

    double a[3] = {1., 5., 3.};
    p->push("A", c2numpy_double(a, 3));
    p->exec("assert (A == array([1., 5., 3.])).all()");

    // integer arrays
    p->exec("a = array(range(20), dtype='int32')");
    p->exec("assert a.strides == (4,)");
    p->exec("b = a[::5]");
    p->exec("assert b.strides == (20,)");
    int *B;
    numpy2c_int_inplace(p->pull("a"), &B, &n);
    _assert(n == 20);
    _assert(
            (B[0] == 0) &&
            (B[1] == 1) &&
            (B[2] == 2) &&
            (B[3] == 3)
           );
    numpy2c_int_inplace(p->pull("b"), &B, &n);
    _assert(n == 4);
    _assert(
            (B[0] == 0) &&
            (B[1] == 5) &&
            (B[2] == 10) &&
            (B[3] == 15)
           );

    int b[3] = {1, 5, 3};
    p->push("B", c2numpy_int(b, 3));
    p->exec("assert (B == array([1, 5, 3])).all()");
}

int main(int argc, char* argv[])
{
    try {
        p = new Python(argc, argv);

        test_basic();
        test_numpy();

        delete p;

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
