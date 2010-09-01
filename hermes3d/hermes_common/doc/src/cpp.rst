.. _cpp_doc:

C++ documentation
=================

This section documents the C++ interface of hermes_common. If you want to use
hermes_common from Python, see the :ref:`Python documentation <python_doc>`.

Below we describe in detail how Hermes Common can be used from C++.

Using Python
------------

Hermes Common provides a ``Python`` C++ class for easy interaction with Python.
Example::

    Python *p = new Python();
    p->push("i", c2py_int(5));
    p->exec("i = i*2");
    int i = py2c_int(p->pull("i"));
    _assert(i == 10);
    delete p;

As you can see, you just instantiate the ``Python`` class, which automatically
initializes Python and sets everything up, then you convert from C++ to Python
using the ``c2py_*`` functions, push it to the Python namespace using
``p->push()``, do whatever you have to do in Python using ``p->exec()``, then
pull it back to C++ using ``p->pull()`` and finally converting from Python to
C++ using the ``py2c_*`` functions.

Python objects are properly allocated and deallocated (e.g. no memory leaks),
so you don't have to worry about reference counting, as long as your code looks
like the example above.

Each ``Python`` instance has it's own namespace, as demonstrated in this
example::

    Python *p1 = new Python();
    Python *p2 = new Python();
    p1->push("i", c2py_int(5));
    p2->push("i", c2py_int(6));
    p1->exec("i = i*2");
    p2->exec("i = i*2");
    int i1 = py2c_int(p1->pull("i"));
    int i2 = py2c_int(p2->pull("i"));
    _assert(i1 == 10);
    _assert(i2 == 12);
    delete p1;
    delete p2;

Matrix Support
--------------

CooMatrix
~~~~~~~~~

Example of creating a real CooMatrix::

    CooMatrix m(5);
    m.add(1, 3, 3.5);
    m.add(2, 3, 4.5);
    m.add(3, 4, 1.5);
    m.add(4, 2, 1.5);
    m.add(2, 3, 1);
    m.print();

and then converting it to CSRMatrix::

    CSRMatrix n1(&m);
    n1.print();

or CSCMatrix::

    CSCMatrix n2(&m);
    n2.print();


Here is how to create a complex CooMatrix and then convert to a complex
CSRMatrix and CSCMatrix::

    CooMatrix m(5, true);
    m.add(1, 3, cplx(3.5));
    m.add(2, 3, cplx(4.5));
    m.add(3, 4, cplx(1.5));
    m.add(4, 2, cplx(1.5));
    m.add(2, 3, cplx(1));
    m.print();

    // convert from COO
    CSRMatrix n1(&m);
    n1.print();
    CSCMatrix n2(&m);
    n2.print();

CSRMatrix and CooMatrix
~~~~~~~~~~~~~~~~~~~~~~~

These matrices can't be constructed from scratch, but you first have to create
some other matrix first (in most cases CooMatrix) and convert it to CSR
afterwards::

    CooMatrix m(5);
    m.add(1, 3, 3.5);
    m.add(2, 3, 4.5);
    m.add(3, 4, 1.5);
    m.add(4, 2, 1.5);
    m.add(2, 3, 1);

    CSRMatrix n1(&m);
    n1.print();

Solvers
-------

Hermes Common provides the following solvers::

    // Dense solvers:
    void solve_linear_system_numpy(Matrix *mat, double *res);
    void solve_linear_system_numpy(Matrix *mat, cplx *res);
    void solve_linear_system_dense_lu(Matrix *mat, double *res);

    // Sparse solvers:
    void solve_linear_system_scipy_umfpack(Matrix *mat, double *res);
    void solve_linear_system_scipy_umfpack(Matrix *mat, cplx *res);
    void solve_linear_system_scipy_cg(Matrix *mat, double *res);
    void solve_linear_system_scipy_gmres(Matrix *mat, double *res);
    int solve_linear_system_cg(Matrix* A, double *x,
                               double matrix_solver_tol,
                               int matrix_solver_maxiter);

They are mostly implemented in SciPy or NumPy (except
``solve_linear_system_dense_lu`` and ``solve_linear_system_cg`` that are
actually implemented in Hermes Common itself in C++) and the implementation
just uses the ``Python`` class to call the corresponding SciPy/NumPy function.
As you can see, all of them accept the abstract Matrix class, so you can supply
a matrix in any format you want and it will be automatically converted (if
needed) to the format that the solver needs (e.g. umfpack needs CSCMatrix,
numpy needs DenseMatrix and so on).

Example::

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

You can use any other solver instead of ``solve_linear_system_dense_lu``.

The solvers handle nonsymmetric matrices just fine, for example::

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

As well as complex numbers, the following example shows that both numpy and
umfpack return the same (complex) solution::

    CooMatrix A(2, true);
    A.add(0, 0, cplx(1, 1));
    A.add(0, 1, cplx(2, 2));
    A.add(1, 0, cplx(3, 3));
    A.add(1, 1, cplx(4, 4));

    cplx res[2];

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
