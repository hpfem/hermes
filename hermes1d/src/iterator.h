// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _ITERATOR_H_
#define _ITERATOR_H_

#include "common.h"
#include "mesh.h"
#include <stack>

class Iterator {
public:
  Iterator(Mesh *mesh) 
  {
    this->mesh = mesh;
    current_coarse_elem_index = -1;
  }
  Mesh *mesh;
  std::stack <Element*>S;  
  void reset();
  Element *first_active_element();
  Element *next_active_element();
  Element *last_active_element();
  int current_coarse_elem_index;
};

#endif
