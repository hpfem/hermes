// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __H1D_WEAKFORM_H
#define __H1D_WEAKFORM_H

/// \brief Represents the weak formulation of a problem.
///
/// The WeakForm class represents the weak formulation of a system of linear PDEs.
/// The number of equations ("neq") in the system is fixed and is passed to the constructor.
/// The weak formulation of the system A(U,V) = L(V) has a block structure. A(U,V) is
/// a (neq x neq) matrix of bilinear forms a_mn(u,v) and L(V) is a neq-component vector
/// of linear forms l(v). U and V are the vectors of basis and test functions.
///
///
///

class HERMES_API WeakForm
{
public:

  WeakForm(int neq = 1, bool mat_free = false);

  // general case
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

  // general case
  void add_matrix_form(int i, int j, matrix_form fn, int marker=ANY);
  void add_matrix_form(matrix_form fn, int marker=ANY); // single equation case
  
  void add_vector_form(int i, vector_form fn, int marker=ANY);
  void add_vector_form(vector_form fn, int marker=ANY); // single equation case
  
  void add_matrix_form_surf(int i, int j, matrix_form_surf fn, int bdy_index);
  void add_matrix_form_surf(matrix_form_surf fn, int bdy_index); // single equation case
  
  void add_vector_form_surf(int i, vector_form_surf fn, int bdy_index);
  void add_vector_form_surf(vector_form_surf fn, int bdy_index); // single equation case

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

  /// Returns the number of equations
  int get_neq() { return neq; }

  /// Internal. Used by DiscreteProblem to detect changes in the weakform.
  int get_seq() const { return seq; }

  bool is_matrix_free() { return is_matfree; }

protected:

  int neq;
  int seq;
  bool is_matfree;
};

#endif
