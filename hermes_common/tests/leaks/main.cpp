#include <iostream>
#include <stdexcept>

#include "python_api.h"
#include "my_api2_api.h"

#include "cpp_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

Python *python;

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

// Test the creation of a big numpy array
void use_big_numpy_array1(int size)
{
    Python *p = new Python();
    p->push("i", c2py_int(size));
    p->exec("import numpy");
    p->exec("A = numpy.zeros((i,), dtype='double')");
    p->exec("print A.shape, 'size (MB):', A.nbytes/1024.**2");
    delete p;
}

void test_leaks1()
{
    // If the matrices are not deallocated, this should take ~19GB:
    for (int i=0; i<5; i++)
        // Roughly 3.8GB:
        use_big_numpy_array1((int)(5*pow(10,8)));
}

class MyCallback: public CppCallback {
    public:
        MyCallback(): CppCallback() {
            this->called = 0;
        }
        virtual void event(const char *msg) {
            printf("MyCallback event(%s) got called.\n", msg);
            this->called = 1;
        }
        int called;
};

void test_dealloc1()
{
    Python *p = new Python();
    if (import_my_api2())
        throw std::runtime_error("my_api failed to import.");
    MyCallback *c = new MyCallback();
    p->push("c", c2py_CppCallback(c));
    _assert(c->called == 0);
    p->exec("from dealloc import A");
    p->exec("a = A(c)");
    p->exec("print 'A() was created'");
    _assert(c->called == 0);
    delete p;
    // This tests that the __del__ method of A() was callled:
    _assert(c->called == 1);
}

void test_dealloc2()
{
    MyCallback *c = new MyCallback();
    {
        Python p = Python();
        p.push("c", c2py_CppCallback(c));
        _assert(c->called == 0);
        p.exec("from dealloc import A");
        p.exec("a = A(c)");
        p.exec("print 'A() was created'");
        _assert(c->called == 0);
    }
    // This tests that the __del__ method of A() was callled:
    _assert(c->called == 1);
}

int main(int argc, char* argv[])
{
    try {
        // Enable this by hand if you want to test with big matrices (might
        // consume all your memory if something goes wrong):
        //test_leaks1();
        test_dealloc1();
        test_dealloc2();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
