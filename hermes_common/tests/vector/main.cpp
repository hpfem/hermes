#include <iostream>
#include <stdexcept>

#include "matrix.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

#define EPS 1e-12

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_vector1()
{
    AVector m(3);
    _assert(m.get_size() == 3);
    m.add(0, 3.5);
    m.add(1, 4.5);
    m.add(2, 1.5);
    m.add(1, 1);
    m.add(0, 2.5);
    m.print();
    _assert(fabs(m.get(0) - 6) < EPS);
    _assert(fabs(m.get(1) - 5.5) < EPS);
    _assert(fabs(m.get(2) - 1.5) < EPS);

    double *v = m.get_c_array();
    _assert(fabs(v[0] - 6) < EPS);
    _assert(fabs(v[1] - 5.5) < EPS);
    _assert(fabs(v[2] - 1.5) < EPS);
}

void test_vector2()
{
    AVector m(5, true);
    m.add(1, cplx(2.3, 3.5));
    m.add(2, cplx(1.2, 4.5));
    m.add(3, cplx(2, 1.5));
    m.add(4, cplx(4.3, 1.5));
    m.add(2, cplx(4.1, 1));
    m.print();
    _assert(fabs(m.get_cplx(0).real() - 0) < EPS);
    _assert(fabs(m.get_cplx(0).imag() - 0) < EPS);
    _assert(fabs(m.get_cplx(1).real() - 2.3) < EPS);
    _assert(fabs(m.get_cplx(1).imag() - 3.5) < EPS);
    _assert(fabs(m.get_cplx(2).real() - 5.3) < EPS);
    _assert(fabs(m.get_cplx(2).imag() - 5.5) < EPS);
    _assert(fabs(m.get_cplx(3).real() - 2) < EPS);
    _assert(fabs(m.get_cplx(3).imag() - 1.5) < EPS);
    _assert(fabs(m.get_cplx(4).real() - 4.3) < EPS);
    _assert(fabs(m.get_cplx(4).imag() - 1.5) < EPS);

    cplx *v = m.get_c_array_cplx();
    _assert(fabs(v[0].real() - 0) < EPS);
    _assert(fabs(v[0].imag() - 0) < EPS);
    _assert(fabs(v[1].real() - 2.3) < EPS);
    _assert(fabs(v[1].imag() - 3.5) < EPS);
    _assert(fabs(v[2].real() - 5.3) < EPS);
    _assert(fabs(v[2].imag() - 5.5) < EPS);
    _assert(fabs(v[3].real() - 2) < EPS);
    _assert(fabs(v[3].imag() - 1.5) < EPS);
    _assert(fabs(v[4].real() - 4.3) < EPS);
    _assert(fabs(v[4].imag() - 1.5) < EPS);
}

#include "python_api.h"

void test_vector3()
{
    AVector m(3);
    _assert(m.get_size() == 3);
    m.add(0, 3.5);
    m.add(1, 4.5);
    m.add(2, 1.5);
    m.add(1, 1);
    m.add(0, 2.5);

    Python p;
    p.push("m", c2py_AVector(&m));
    p.exec("print m");
    p.exec("d = m.to_numpy()");
    p.exec("print d");
    p.exec("eps = 1e-10");
    p.exec("assert abs(d[0]-6) < eps");
    p.exec("assert abs(d[1]-5.5) < eps");
    p.exec("assert abs(d[2]-1.5) < eps");
}

int main(int argc, char* argv[])
{
    try {
        test_vector1();
        test_vector2();
        test_vector3();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
