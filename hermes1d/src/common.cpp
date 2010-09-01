// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "common.h"

void error(const char *msg)
{
  printf("Error: %s\n", msg);
  exit(0);
};

void error(const char *msg, const char *msg1)
{
  printf("Error: %s %s\n", msg, msg1);
  exit(0);
};

void info(const char *msg, const char *msg1)
{
  printf("Info: %s %s\n", msg, msg1);
};

void warning(const char *msg)
{
  printf("Warning: %s\n", msg);
};

void intro() {
  printf("\n-------------------------------------------\n");
  printf("   This is Hermes1D - a free ODE solver\n");
  printf(" based on the hp-FEM and Newton's method,\n");
  printf("   developed by the hp-FEM group at UNR\n");
  printf("  and distributed under the BSD license.\n");
  printf(" For more details visit http://hpfem.org/.\n");
  printf("-------------------------------------------\n");
}
