#include "stdcython.h"

PyObject* global_empty_tuple;

void init_global_empty_tuple() {
  _CALLED_ONLY_ONCE;

  global_empty_tuple = PyTuple_New(0);
}
