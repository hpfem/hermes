#include <iostream>
#include <stdexcept>

#include "python_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

Python *python;

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

// Test the global python instance
void test_basic1()
{
    python->push("i", c2py_int(5));
    python->exec("i = i*2");
    int i = py2c_int(python->pull("i"));
    _assert(i == 10);
}

// Test the local Python() instance
void test_basic2()
{
    Python *p = new Python();
    p->push("i", c2py_int(5));
    p->exec("i = i*2");
    int i = py2c_int(p->pull("i"));
    _assert(i == 10);
    delete p;
}

// Test initialization/destruction of two Python() instances
void test_basic3()
{
    Python *p1 = new Python();
    p1->push("i", c2py_int(5));
    p1->exec("i = i*2");
    int i = py2c_int(p1->pull("i"));
    _assert(i == 10);
    delete p1;

    Python *p2 = new Python();
    p2->push("i", c2py_int(5));
    p2->exec("i = i*2");
    i = py2c_int(p2->pull("i"));
    _assert(i == 10);
    delete p2;
}

// Test that each Python() instance has it's own namespace
void test_basic4()
{
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
}

int main(int argc, char* argv[])
{
    try {
        python = new Python(argc, argv);

        test_basic1();
        test_basic2();
        test_basic3();
        test_basic4();

        delete python;

        // this tests that all Python() instances were uninitialized properly
        test_basic3();
        test_basic4();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
