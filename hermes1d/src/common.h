// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES1D_COMMON_H
#define __HERMES1D_COMMON_H

#include <stdexcept>
#include <stdlib.h>
#include <stdio.h>

// print general debug information
#define DEBUG 0

// printf debug information about the stiffness/Jacobi matrix
#define DEBUG_MATRIX 0

#define BOUNDARY_LEFT 0
#define BOUNDARY_RIGHT 1

// for material flags
const int ANY = -1234;

// Up to 100 is currently implemented.
// When you change this, run tests.
const int MAX_P = 30;                  // max poly degree allowed in elements
                                       // WARNING: projections taking place in 
                                       // transfer_solution()
// Up to 200 is currently implemented:
// When you change this, run tests.
const int MAX_QUAD_ORDER = 200;        // max order of Gaussian quadrature implemented
const int MAX_QUAD_PTS_NUM = 101;      // max number of quadrature points

const int MAX_CAND_NUM = 100;          // maximum allowed number of hp-refinement
                                       // candidates of an element

const int MAX_ELEM_NUM = 10000;        // maximum number of elements
const int MAX_N_DOF = 10000;           // maximum number of degrees of freedom
const int MAX_EQN_NUM = 10;            // maximum number of equations in the system
const int MAX_SLN_NUM = 5;             // maximum number of solutions (not to be confused 
                                       // components - every solution can have multiple 
                                       // components)
const int MAX_PLOT_PTS_NUM = 501;      // finest plotting subdivision 
const int MAX_COEFFS_NUM = MAX_P + 1;  // this is the maximum number of 
                                       // polynomial coefficients
const int MAX_STRING_LENGTH = 100;     // maximum string length 

typedef double (*exact_sol_type)(double x, 
        double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]);

void error(const char *msg);
void error(const char *msg, const char *msg1);
void info(const char *msg, const char *msg1);
void warning(const char *msg);

// types
typedef double scalar;
typedef double double2[2];
typedef int int2[2];
typedef int int3[3];
typedef double (*shape_fn_t)(double);

// auxiliary functions
void intro();
#define MEM_CHECK(var) if (var == NULL) { printf("Out of memory."); exit(1); }
#define verbose(msg)
#define warn(msg)

#endif
