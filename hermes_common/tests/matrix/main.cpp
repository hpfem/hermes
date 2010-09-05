#include <iostream>
#include <stdexcept>

#include "matrix.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_matrix1()
{
    CooMatrix m(3);
    m.add(0, 2, 3.5);
    m.add(1, 2, 4.5);
    m.add(2, 1, 1.5);
    m.add(1, 2, 1);
    m.add(0, 0, 2.5);
    m.print();

    printf("----\n");
    // convert from COO
    CSRMatrix n1(&m);
    n1.print();
    CSCMatrix n2(&m);
    n2.print();

    // convert CSR <-> CSC
    CSRMatrix n3(&n2);
    n3.print();
    CSCMatrix n4(&n1);
    n4.print();
}

void test_matrix2()
{
    CooMatrix m(5);
    m.add(1, 3, 3.5);
    m.add(2, 3, 4.5);
    m.add(3, 4, 1.5);
    m.add(4, 2, 1.5);
    m.add(2, 3, 1);
    m.print();

    Matrix *_m = &m;

    printf("----\n");
    // convert from COO
    CSRMatrix n1(_m);
    n1.print();    
    CSCMatrix n2(_m);
    n2.print();

    _m = &n2;
    // convert CSR <-> CSC
    CSRMatrix n3(_m);
    n3.print();
    _m = &n1;
    CSCMatrix n4(_m);
    n4.print();
}

void test_matrix3()
{
    CooMatrix m(5, true);
    m.add(1, 3, cplx(2.3, 3.5));
    m.add(2, 3, cplx(1.2, 4.5));
    m.add(3, 4, cplx(2, 1.5));
    m.add(4, 2, cplx(4.3, 1.5));
    m.add(2, 3, cplx(4.1, 1));
    m.print();

    printf("----\n");
    // convert from COO
    CSRMatrix n1(&m);
    n1.print();
    CSCMatrix n2(&m);
    n2.print();

    // convert CSR <-> CSC
    CSRMatrix n3(&n2);
    n3.print();
    CSCMatrix n4(&n1);
    n4.print();

}

void test_matrix4()
{
    DenseMatrix m(3);
    m.add(0, 2, 3.5);
    m.add(1, 2, 4.5);
    m.add(2, 1, 1.5);
    m.add(1, 2, 1);
    m.add(0, 0, 2.5);
    printf("Dense Matrix\n");
    m.print();

    printf("----\n");
    // convert to CSR and CSC
    CSRMatrix n1(&m);
    n1.print();
    /*
    FIXME: DenseMatrix -> CSCMatrix segfaults...:
    CSCMatrix n2(&m);
    n2.print();
    */
}

#include "python_api.h"

void test_matrix5()
{
    DenseMatrix m(3);
    m.add(0, 2, 3.5);
    m.add(1, 2, 4.5);
    m.add(2, 1, 1.5);
    m.add(1, 2, 1);
    m.add(0, 0, 2.5);
    printf("Dense Matrix\n");

    Python p;
    p.push("m", c2py_DenseMatrix(&m));
    p.exec("print m");
    p.exec("d = m.to_numpy()");
    p.exec("print d");
    p.exec("eps = 1e-10");
    p.exec("assert abs(d[0, 0]-2.5) < eps");
    p.exec("assert abs(d[0, 2]-3.5) < eps");
    p.exec("assert abs(d[1, 2]-5.5) < eps");
    p.exec("assert abs(d[2, 1]-1.5) < eps");
    p.exec("assert abs(d[0, 1]-0.0) < eps");
}

int main(int argc, char* argv[])
{
    try {
        test_matrix1();
        test_matrix2();
        test_matrix3();
        test_matrix4();
        test_matrix5();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
