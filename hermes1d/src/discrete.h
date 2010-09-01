// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _DISCRETE_H_
#define _DISCRETE_H_

#include <vector>

#include "mesh.h"
#include "quad_std.h"
#include "legendre.h"
#include "lobatto.h"
#include "matrix.h"
#include "iterator.h"

typedef double (*matrix_form) (int num, double *x, double *weights,
        double *u, double *dudx, double *v, double *dvdx, 
        double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], void *user_data);

typedef double (*vector_form) (int num, double *x, double *weights,
        double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
               double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
               double *v, double *dvdx,
        void *user_data);

typedef double (*matrix_form_surf) (double x, double u, double dudx, 
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data);

typedef double (*vector_form_surf) (double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v, double dvdx,
        void *user_data);

class DiscreteProblem {

public:
    DiscreteProblem();
    void add_matrix_form(int i, int j, matrix_form fn, int marker=ANY);
    void add_vector_form(int i, vector_form fn, int marker=ANY);
    void add_matrix_form_surf(int i, int j, matrix_form_surf fn, int bdy_index);
    void add_vector_form_surf(int i, vector_form_surf fn, int bdy_index);
    // c is solution component
    void process_vol_forms(Mesh *mesh, Matrix *mat, double *res, 
                           int matrix_flag);
    // c is solution component
    void process_surf_forms(Mesh *mesh, Matrix *mat, double *res, 
                            int matrix_flag, int bdy_index);
    void assemble(Mesh *mesh, Matrix *mat, double *res, int matrix_flag);
    void assemble_matrix_and_vector(Mesh *mesh, Matrix *mat, double *res); 
    void assemble_matrix(Mesh *mesh, Matrix *mat);
    void assemble_vector(Mesh *mesh, double *res);

private:
	struct MatrixFormVol {
		int i, j;
		matrix_form fn;
	        int marker;
	};
	struct MatrixFormSurf {
		int i, j, bdy_index;
		matrix_form_surf fn;
	};
	struct VectorFormVol {
		int i;
		vector_form fn;
	        int marker;
	};
	struct VectorFormSurf {
		int i, bdy_index;
		vector_form_surf fn;
	};
	std::vector<MatrixFormVol> matrix_forms_vol;
	std::vector<MatrixFormSurf> matrix_forms_surf;
	std::vector<VectorFormVol> vector_forms_vol;
	std::vector<VectorFormSurf> vector_forms_surf;
};

// return coefficients for all shape functions on the element m,
// for all solution components
void calculate_elem_coeffs(Element *e, double **coeffs, int n_eq);

void element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int *num);

void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der);

void element_shapefn_point(double x_ref, double a, double b, 
			   int k, double &val, double &der);

void newton(DiscreteProblem *dp, Mesh *mesh, 
            CommonSolver *solver,
            double newton_tol, int newton_maxiter,
            bool verbose=true);

void jfnk_cg(DiscreteProblem *dp, Mesh *mesh,
             double matrix_solver_tol, int matrix_solver_maxiter,  
	     double jfnk_epsilon, double jfnk_tol, int jfnk_maxiter, bool verbose=true);



#endif
