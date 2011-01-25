#include <stdexcept>
#include <csignal>

#include "python_api.h"
#include "python_engine_api.h"

static int python_count=0;

Python::Python()
{
    this->_init(-1, NULL);
}

Python::Python(int argc, char* argv[])
{
    this->_init(argc, argv);
}

PyMODINIT_FUNC initpython_engine(void); /*proto*/

void exit_program(int sig) {
    printf("\nCTRL-C received, aborting...\n");
    abort();
}

void Python::_init(int argc, char* argv[])
{
    python_count++;
    if (python_count == 1) {
        Py_Initialize();
        if (argc >= 0)
            PySys_SetArgv(argc, argv);
        // We need to install our own signal, because Python redefines it and
        // then the C++ program can't be killed with CTRL-C:
        signal(SIGINT, exit_program);
        // This initializes the Python module using Python C / API:
        initpython_engine();
        // This imports the module using Cython functionality, so that we can
        // use the Cython api functions
        if (import_python_engine())
            throw std::runtime_error("python_engine failed to import.");
    }
    this->_namespace = namespace_create();
}

Python::~Python()
{
    // Free the namespace. This frees all the dictionary items, so if there
    // are some numpy arrays (or your classes) in the namespace, they will be
    // deallocated at this time.
    Py_DECREF(this->_namespace);

    // The code below would free the interpreter if this was the last instance
    // using it. However, it is currently disabled, because the numpy package
    // segfaults when imported again; also the PYTHONPATH is set only once if
    // python_count is never decreased (which is what we want).
    /*python_count--;
    if (python_count == 0) {
        Py_Finalize();
    }
    */
}

void Python::print_namespace()
{
    namespace_print(_namespace);
}

void Python::exec(const std::string &text)
{
    run_cmd(text.c_str(), this->_namespace);
}

void Python::push(const std::string &name, PyObject *o)
{
    namespace_push(this->_namespace, name.c_str(), o);
    // namespace_push() is a regular Cython function and
    // as such, it increfs the object "o" before storing it in the namespace,
    // but we want to steal the reference, so we decref it here (there is still
    // at least one reference stored in the dictionary this->_namespace, so
    // it's safe). This is so that
    //     this->push("i", c2py_int(5));
    // doesn't leak (c2py_int() creates a python reference and push() destroys
    // this python reference)
    Py_DECREF(o);
}

PyObject *Python::pull(const std::string &name)
{
    PyObject *tmp = namespace_pull(this->_namespace, name.c_str());
    // namespace_pull() is a regular Cython function and
    // as such, it increfs the result before returning it, but we only want to
    // borrow a reference, so we decref it here (there is still at least one
    // reference stored in the dictionary this->_namespace, so it's safe)
    // This is so that
    //     int i = py2c_int(this->pull("i"));
    // doesn't leak (pull() borrows the reference, py2c_int() doesn't do
    // anything with the reference, so no leak nor segfault happens)
    Py_DECREF(tmp);
    return tmp;
}

void Python::push_int(const std::string &name, int i)
{
    this->push(name, c2py_int(i));
}

int Python::pull_int(const std::string &name)
{
    return py2c_int(this->pull(name));
}

void Python::push_double(const std::string &name, double i)
{
    this->push(name, c2py_double(i));
}

double Python::pull_double(const std::string &name)
{
    return py2c_double(this->pull(name));
}

void Python::push_numpy_double(const std::string &name, double *A, int n)
{
    this->push(name, c2numpy_double(A, n));
}

void Python::push_numpy_double_inplace(const std::string &name,
        double *A, int n)
{
    this->push(name, c2numpy_double_inplace(A, n));
}

void Python::pull_numpy_double_inplace(const std::string &name,
        double **A, int *n)
{
    numpy2c_double_inplace(this->pull(name), A, n);
}

void Python::push_numpy_int(const std::string &name, int *A, int n)
{
    this->push(name, c2numpy_int(A, n));
}

void Python::push_numpy_int_inplace(const std::string &name, int *A, int n)
{
    this->push(name, c2numpy_int_inplace(A, n));
}

void Python::pull_numpy_int_inplace(const std::string &name,
        int **A, int *n)
{
    numpy2c_int_inplace(this->pull(name), A, n);
}

void Python::push_str(const std::string &name, const std::string &s)
{
    this->push(name, c2py_str(s.c_str()));
}

std::string Python::pull_str(const std::string &name)
{
    return py2c_str(this->pull(name));
}
