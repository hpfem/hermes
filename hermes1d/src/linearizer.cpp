// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "mesh.h"
#include "linearizer.h"
#include "iterator.h"

// Evaluate (vector-valued) approximate solution at reference 
// point 'x_ref' in element 'm'. Here 'y' is the global vector 
// of coefficients. The result is a vector of length mesh->n_eq
void Linearizer::eval_approx(int sln, Element *e, double x_ref,
                             double *x_phys, double *val) 
{
  int n_eq = this->mesh->get_n_eq();
  for(int c=0; c<n_eq; c++) { // loop over solution components
    val[c] = 0;
    for(int i=0; i <= e->p; i++) { // loop over shape functions
      if(e->dof[c][i] >= 0) val[c] += 
                  e->coeffs[sln][c][i]*lobatto_val_ref(x_ref, i);
    }
  }
  double a = e->x1;
  double b = e->x2;
  *x_phys = (a+b)/2 + x_ref*(b-a)/2;
  return;
}
// default is for sln=0
void Linearizer::eval_approx(Element *e, double x_ref,
                             double *x_phys, double *val) 
{
  eval_approx(0, e, x_ref, x_phys, val); 
}

// Plot solution in Gnuplot format
void Linearizer::plot_solution(const char *out_filename, 
                               int plotting_elem_subdivision)
{
    int n_eq = this->mesh->get_n_eq();
    FILE *f[MAX_EQN_NUM];
    char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
    for(int c=0; c<n_eq; c++) {
        if(n_eq == 1)
            sprintf(final_filename[c], "%s", out_filename);
        else
            sprintf(final_filename[c], "%s_%d", out_filename, c);
        f[c] = fopen(final_filename[c], "wb");
        if(f[c] == NULL) error("problem opening file in plot_solution().");
        int n;
        double *x, *y;
        this->get_xy_mesh(c, plotting_elem_subdivision, &x, &y, &n);
        for (int i=0; i < n; i++)
            fprintf(f[c], "%g %g\n", x[i], y[i]);
        fprintf(f[c], "\n");
        delete[] x;
        delete[] y;
        printf("Output written to %s.\n", final_filename[c]);
        fclose(f[c]);
    }
}

// Plot solution stored in the ref_elem_pairs[] array 
// in Gnuplot format
void Linearizer::plot_ref_elem_pairs(ElemPtr2* ref_elem_pairs, 
                                     const char *out_filename, 
                                     int plotting_elem_subdivision)
{
    int n_eq = this->mesh->get_n_eq();
    FILE *f[MAX_EQN_NUM];
    char final_filename[MAX_EQN_NUM][MAX_STRING_LENGTH];
    for(int c=0; c<n_eq; c++) {
        if(n_eq == 1)
            sprintf(final_filename[c], "%s", out_filename);
        else
            sprintf(final_filename[c], "%s_%d", out_filename, c);
        f[c] = fopen(final_filename[c], "wb");
        if(f[c] == NULL) error("problem opening file in plot_solution().");
        int n;
        double *x, *y;
        this->get_xy_ref_array(c, ref_elem_pairs, plotting_elem_subdivision, &x, &y, &n);
        for (int i=0; i < n; i++)
            fprintf(f[c], "%g %g\n", x[i], y[i]);
        fprintf(f[c], "\n");
        delete[] x;
        delete[] y;
        printf("Output written to %s.\n", final_filename[c]);
        fclose(f[c]);
    }
}

// Plot solution in Gnuplot format, as a trajectory where solution[comp_x]
// is used as x-coordinate and solution[comp_y] as y-coordinate
void Linearizer::plot_trajectory(FILE *f,
			         int comp_x, int comp_y, 
                                 int plotting_elem_subdivision)
{
    int n1, n2;
    double *x1, *y1, *x2, *y2;
    this->get_xy_mesh(comp_x, plotting_elem_subdivision, &x1, &y1, &n1);
    this->get_xy_mesh(comp_y, plotting_elem_subdivision, &x2, &y2, &n2);
    if (n1 != n2) error("internal: n1 != n2 in plot_trajectory().");
    for (int i=0; i < n1; i++) fprintf(f, "%g %g\n", y1[i], y2[i]);
    fprintf(f, "\n");
    delete[] x1;
    delete[] y1;
    delete[] x2;
    delete[] y2;
}


// Returns pointers to x and y coordinates in **x and **y
// you should free it yourself when you don't need it anymore
// comp --- which component you want to process
// plotting_elem_subdivision --- the number of subdivision of the element
// x, y --- the doubles list of x,y
// n --- the number of points

void Linearizer::get_xy_mesh(int comp,
                        int plotting_elem_subdivision,
                        double **x, double **y, int *n)
{
    int n_eq = this->mesh->get_n_eq();
    int n_active_elem = this->mesh->get_n_active_elem();
    Iterator *I = new Iterator(this->mesh);

    *n = n_active_elem * (plotting_elem_subdivision+1);
    double *x_out = new double[*n];
    double *y_out = new double[*n];

    // FIXME:
    if(n_eq > MAX_EQN_NUM) {
      printf("n_eq = %d\n", n_eq);
        error("number of equations too high in plot_solution().");
    }
    // FIXME
    if(plotting_elem_subdivision > MAX_PLOT_PTS_NUM)
        error("plotting_elem_subdivision too high in plot_solution().");
    double phys_u_prev[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
    double phys_du_prevdx[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
        
    Element *e;
    int counter = 0;
    while ((e = I->next_active_element()) != NULL) {
        //printf("linearizer: in element (%g, %g)\n", e->x1, e->x2);
        if (counter >= n_active_elem) {
	  printf("n_active_elem = %d\n", n_active_elem);
	  printf("counter = %d\n", counter);
            error("Internal error: wrong n_active_elem");
        }

        double x_phys[MAX_PLOT_PTS_NUM];
        double h = (e->x2 - e->x1)/plotting_elem_subdivision;

        for (int j=0; j<plotting_elem_subdivision+1; j++)
            x_phys[j] = e->x1 + j*h;
        e->get_solution_plot(x_phys, plotting_elem_subdivision+1, 
                phys_u_prev, phys_du_prevdx);
        double a = e->x1;
        double b = e->x2;
        for (int j=0; j<plotting_elem_subdivision+1; j++) {
            x_out[counter*(plotting_elem_subdivision+1) + j] = x_phys[j];
            y_out[counter*(plotting_elem_subdivision+1) + j] =
                phys_u_prev[comp][j];
        }
        counter++;
    }
    *x = x_out;
    *y = y_out;
    delete I;
}

void Linearizer::get_xy_ref_array(int comp, ElemPtr2* ref_elem_pairs,
                        int plotting_elem_subdivision,
                        double **x, double **y, int *n)
{
    int n_eq = this->mesh->get_n_eq();
    int n_active_elem = this->mesh->get_n_active_elem();

    *n = 2 * n_active_elem * (plotting_elem_subdivision+1);
    double *x_out = new double[*n];
    double *y_out = new double[*n];

    // FIXME:
    if(n_eq > MAX_EQN_NUM) {
      printf("n_eq = %d\n", n_eq);
        error("number of equations too high in plot_solution().");
    }
    // FIXME
    if(plotting_elem_subdivision > MAX_PLOT_PTS_NUM)
        error("plotting_elem_subdivision too high in plot_solution().");
    double phys_u_prev[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];
    double phys_du_prevdx[MAX_EQN_NUM][MAX_PLOT_PTS_NUM];

    for (int i=0; i < n_active_elem; i++) {
      for (int m=0; m < 2; m++) {
        Element *e = ref_elem_pairs[i][m];
        double x_phys[MAX_PLOT_PTS_NUM];
        double h = (e->x2 - e->x1)/plotting_elem_subdivision;

        for (int j=0; j<plotting_elem_subdivision+1; j++)
            x_phys[j] = e->x1 + j*h;
        e->get_solution_plot(x_phys, plotting_elem_subdivision+1, 
                phys_u_prev, phys_du_prevdx);
        double a = e->x1;
        double b = e->x2;
        for (int j=0; j<plotting_elem_subdivision+1; j++) {
	  x_out[(2*i+m)*(plotting_elem_subdivision + 1) + j] = x_phys[j];
          y_out[(2*i+m)*(plotting_elem_subdivision + 1) + j] =
                phys_u_prev[comp][j];
        }
      }
    }
    *x = x_out;
    *y = y_out;
}


