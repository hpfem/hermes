// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file python_engine_api.h
    \brief Python API engine.
*/
#ifndef __PYX_HAVE_API__python_engine
#define __PYX_HAVE_API__python_engine
#ifdef WITH_PYTHON

#include "Python.h"

static PyObject *(*namespace_create)(void);
static void (*namespace_push)(PyObject *, const char*, PyObject *);
static void (*namespace_print)(PyObject *);
static PyObject *(*namespace_pull)(PyObject *, const char*);
static PyObject *(*c2py_int)(int);
static int (*py2c_int)(PyObject *);
static PyObject *(*c2py_double)(double);
static double (*py2c_double)(PyObject *);
static PyObject *(*c2py_str)(const char*);
static char *(*py2c_str)(PyObject *);
static PyObject *(*c2numpy_int)(int *, int);
static PyObject *(*c2numpy_int_inplace)(int *, int);
static PyObject *(*c2numpy_double)(double *, int);
static PyObject *(*c2numpy_double_inplace)(double *, int);
static void (*numpy2c_int_inplace)(PyObject *, int **, int *);
static void (*numpy2c_double_inplace)(PyObject *, double **, int *);
static void (*run_cmd)(const char*, PyObject *);

#ifndef __PYX_HAVE_API_FUNC_import_module
#define __PYX_HAVE_API_FUNC_import_module

#ifndef __PYX_HAVE_RT_ImportModule
#define __PYX_HAVE_RT_ImportModule
static PyObject *__Pyx_ImportModule(const char *name) {
    PyObject *py_name = 0;
    PyObject *py_module = 0;

    #if PY_MAJOR_VERSION < 3
    py_name = PyString_FromString(name);
    #else
    py_name = PyUnicode_FromString(name);
    #endif
    if (!py_name)
        goto bad;
    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);
    return py_module;
bad:
    Py_XDECREF(py_name);
    return 0;
}
#endif

#endif


#ifndef __PYX_HAVE_RT_ImportFunction
#define __PYX_HAVE_RT_ImportFunction
static int __Pyx_ImportFunction(PyObject *module, const char *funcname, void (**f)(void), const char *sig) {
    PyObject *d = 0;
    PyObject *cobj = 0;
    union {
        void (*fp)(void);
        void *p;
    } tmp;

    d = PyObject_GetAttrString(module, (char *)"__pyx_capi__");
    if (!d)
        goto bad;
    cobj = PyDict_GetItemString(d, funcname);
    if (!cobj) {
        PyErr_Format(PyExc_ImportError,
            "%s does not export expected C function %s",
                PyModule_GetName(module), funcname);
        goto bad;
    }
#if PY_VERSION_HEX >= 0x02070000 && !(PY_MAJOR_VERSION==3&&PY_MINOR_VERSION==0)
    if (!PyCapsule_IsValid(cobj, sig)) {
        PyErr_Format(PyExc_TypeError,
            "C function %s.%s has wrong signature (expected %s, got %s)",
             PyModule_GetName(module), funcname, sig, PyCapsule_GetName(cobj));
        goto bad;
    }
    tmp.p = PyCapsule_GetPointer(cobj, sig);
#else
    {const char *desc, *s1, *s2;
    desc = (const char *)PyCObject_GetDesc(cobj);
    if (!desc)
        goto bad;
    s1 = desc; s2 = sig;
    while (*s1 != '\0' && *s1 == *s2) { s1++; s2++; }
    if (*s1 != *s2) {
        PyErr_Format(PyExc_TypeError,
            "C function %s.%s has wrong signature (expected %s, got %s)",
             PyModule_GetName(module), funcname, sig, desc);
        goto bad;
    }
    tmp.p = PyCObject_AsVoidPtr(cobj);}
#endif
    *f = tmp.fp;
    if (!(*f))
        goto bad;
    Py_DECREF(d);
    return 0;
bad:
    Py_XDECREF(d);
    return -1;
}
#endif

static int import_python_engine(void) {
  PyObject *module = 0;
  module = __Pyx_ImportModule("python_engine");
  if (!module) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_create", (void (**)(void))&namespace_create, "PyObject *(void)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_push", (void (**)(void))&namespace_push, "void (PyObject *, const char*, PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_print", (void (**)(void))&namespace_print, "void (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "namespace_pull", (void (**)(void))&namespace_pull, "PyObject *(PyObject *, const char*)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_int", (void (**)(void))&c2py_int, "PyObject *(int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_int", (void (**)(void))&py2c_int, "int (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_double", (void (**)(void))&c2py_double, "PyObject *(double)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_double", (void (**)(void))&py2c_double, "double (PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2py_str", (void (**)(void))&c2py_str, "PyObject *(const char*)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "py2c_str", (void (**)(void))&py2c_str, "char *(PyObject *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_int", (void (**)(void))&c2numpy_int, "PyObject *(int *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_int_inplace", (void (**)(void))&c2numpy_int_inplace, "PyObject *(int *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_double", (void (**)(void))&c2numpy_double, "PyObject *(double *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "c2numpy_double_inplace", (void (**)(void))&c2numpy_double_inplace, "PyObject *(double *, int)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "numpy2c_int_inplace", (void (**)(void))&numpy2c_int_inplace, "void (PyObject *, int **, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "numpy2c_double_inplace", (void (**)(void))&numpy2c_double_inplace, "void (PyObject *, double **, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "run_cmd", (void (**)(void))&run_cmd, "void (const char*, PyObject *)") < 0) goto bad;
  Py_DECREF(module); module = 0;
  return 0;
  bad:
  Py_XDECREF(module);
  return -1;
}
#endif
#endif
