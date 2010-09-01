#ifndef STDCYTHON_H
#define STDCYTHON_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * A global empty python tuple object. This is used to speed up some
 * python API calls where we want to avoid constructing a tuple every
 * time.
 */

extern PyObject* global_empty_tuple;
void init_global_empty_tuple();

/** Constructs a new object of type zzz_type by calling tp_new
 *  directly, with no arguments.
 */

#define PY_NEW(zzz_type) \
    (((PyTypeObject*)(zzz_type))->tp_new((PyTypeObject*)(zzz_type), global_empty_tuple, NULL))

/**
 * a handy macro to be placed at the top of a function definition
 * below the variable declarations to ensure a function is called once
 * at maximum.
 */
#define _CALLED_ONLY_ONCE static int ncalls = 0; if (ncalls>0) return; else ncalls++

#ifdef __cplusplus
}
#endif

#endif

