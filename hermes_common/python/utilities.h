// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_UTILITIES_H
#define __HERMES_COMMON_UTILITIES_H

#include <stdexcept>

#ifdef __cplusplus
extern "C" {
#endif

void throw_exception(char *text)
{
    throw std::runtime_error(text);
}

#ifdef __cplusplus
}
#endif

#endif
