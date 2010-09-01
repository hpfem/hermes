// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _LINEARIZER_H_
#define _LINEARIZER_H_

#include "common.h"
#include "legendre.h"
#include "lobatto.h"

class Linearizer {
    public:
        Linearizer(Mesh *mesh) {
            this->mesh = mesh;
        }

        // evaluate approximate solution at element 'm' at reference
        // point 'x_ref'. Here 'y' is the global vector of coefficients
        void eval_approx(int sln, Element *e, double x_ref, double *x_phys,
			 double *val);
        void eval_approx(Element *e, double x_ref, double *x_phys, // default for sln=0
			 double *val);

        void plot_solution(const char *out_filename, 
                           int plotting_elem_subdivision=50);
        void plot_ref_elem_pairs(ElemPtr2 *elem_ref_pairs,
                           const char *out_filename,
                           int plotting_elem_subdivision=50);
    
        // plotting trajectory where solution[comp_x] is used as
        // x-coordinate and solution[comp_y] as y-coordinate
        // FIXME: code needs to be fixed to allow
        // plotting_elem_subdivision to be 100 and more
        void plot_trajectory(FILE *f, int comp_x, int comp_y, 
                             int plotting_elem_subdivision=50);

        void get_xy_mesh(int comp, int plotting_elem_subdivision,
                         double **x, double **y, int *n);
        void get_xy_ref_array(int comp, ElemPtr2* elem_ref_pairs, 
                              int plotting_elem_subdivision,
                              double **x, double **y, int *n);

    private:
        Mesh *mesh;
};

#endif
