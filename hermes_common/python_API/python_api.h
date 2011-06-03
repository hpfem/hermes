// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_PYTHON_API_H
#define __HERMES_COMMON_PYTHON_API_H

#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE
#undef HAVE_SYS_TIME_H
#include "Python.h"
#include "compat.h"

/*
    This is a nice C++ Python API and the only header file that you should
    include in your code.

    Use the push_int, pull_int (and similar) methods and then you don't need to
    worry about anything. Here is an example how to use it:

        Python *p = new Python();
        p->push_int("i", 5);
        p->exec("i = i*2");
        int i = p->pull_int("i");
        _assert(i == 10);
        delete p;

    All memory allocation/deallocation as well as Python initialization is
    handled automatically.

*/

class HERMES_API Python {
public:
    Python();
    Python(int argc, char* argv[]);
    ~Python();
    void print_namespace();
    void exec(const std::string &text);

    void push_int(const std::string &name, int i);
    int pull_int(const std::string &name);
    void push_double(const std::string &name, double i);
    double pull_double(const std::string &name);
    void push_str(const std::string &name, const std::string &s);
    std::string pull_str(const std::string &name);
    void push_numpy_double(const std::string &name, double *A, int n);
    void push_numpy_double_inplace(const std::string &name, double *A, int n);
    void pull_numpy_double_inplace(const std::string &name, double **A, int *n);
    void push_numpy_int(const std::string &name, int *A, int n);
    void push_numpy_int_inplace(const std::string &name, int *A, int n);
    void pull_numpy_int_inplace(const std::string &name, int **A, int *n);



    /*
    The following two methods are for advanced users only. Make sure you
    understand the Python reference counting and who is responsible for what
    in terms of memory management before you use them.

    Here is an example how to use them:

        Python *p = new Python();
        p->push("i", c2py_int(5));
        p->exec("i = i*2");
        int i = py2c_int(p->pull("i"));
        _assert(i == 10);
        delete p;

    All memory allocation/deallocation as well as Python initialization is
    handled automatically, you don't have to worry about anything, as long as
    you use it as shown above. See also the implementation of pull/push() for
    more info and comments.

    */

    // pushes the object to the namespace, stealing the reference
    void push(const std::string &name, PyObject *o);
    // pulls the object from the namespace, borrowing the reference
    PyObject *pull(const std::string &name);

private:
    PyObject *_namespace;
    void _init(int argc, char* argv[]);
};

#endif
