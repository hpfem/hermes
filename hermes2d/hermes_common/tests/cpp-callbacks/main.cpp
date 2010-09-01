#include <iostream>
#include <stdexcept>

#include "python_api.h"
#include "my_api_api.h"

#include "cpp_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

Python *python;

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_api1()
{
    Python *p = new Python();
    p->exec("from my_api import CppCallback");
    p->exec("c = CppCallback()");
    p->exec("c.event('test')");
    delete p;
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

void test_api2()
{
    Python *p = new Python();
    if (import_my_api())
        throw std::runtime_error("my_api failed to import.");
    MyCallback *c = new MyCallback();
    p->push("c", c2py_CppCallback(c));
    _assert(c->called == 0);
    // This calls MyCallback.event():
    p->exec("c.event('test')");
    _assert(c->called == 1);
    delete p;
    //delete c;
}

void test_callback()
{
    Python *p = new Python();
    p->exec("from callback import A");
    p->exec("a = A()");
    p->exec("print 'A() was created'");
    delete p;
}

int main(int argc, char* argv[])
{
    try {
        test_api1();
        test_api2();
        test_callback();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
