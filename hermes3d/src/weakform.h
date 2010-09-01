// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// This file was taken from hermes2d and adjusted for hermes3d
//

#ifndef _WEAKFORM_H_
#define _WEAKFORM_H_

#include "function.h"
#include "forms.h"
#include "tuple.h"

// Bilinear form symmetry flag, see WeakForm::add_matrix_form
enum SymFlag {
	ANTISYM = -1,
	UNSYM = 0,
	SYM = 1
};

/// Matrix and vector forms.
typedef scalar (*matrix_form_val_t)(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<double> *vi,
	        fn_t<double> *vj, geom_t<double> *e, user_data_t<scalar> *);
typedef ord_t (*matrix_form_ord_t)(int n, double *wt, fn_t<ord_t> *u_ext[], fn_t<ord_t> *vi,
	       fn_t<ord_t> *vj, geom_t<ord_t> *e, user_data_t<ord_t> *);
typedef scalar (*vector_form_val_t)(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<double> *vi,
	        geom_t<double> *e, user_data_t<scalar> *);
typedef ord_t (*vector_form_ord_t)(int n, double *wt, fn_t<ord_t> *u_ext[], fn_t<ord_t> *vi,
	       geom_t<ord_t> *e, user_data_t<ord_t> *);

/// Represents the weak formulation of a problem.
///
/// The WeakForm class represents the weak formulation of a system of linear PDEs.
/// The number of equations ("neq") in the system is fixed and is passed to the constructor.
/// The weak formulation of the system A(U,V) = L(V) has a block structure. A(U,V) is
/// a (neq x neq) matrix of bilinear forms a_mn(u,v) and L(V) is a neq-component vector
/// of linear forms l(v). U and V are the vectors of basis and test functions.
///
///
/// @ingroup assembling
class WeakForm {

public:
	WeakForm(int neq, bool mat_free = false);
	WeakForm(bool mat_free = false);            // single equation case
	virtual ~WeakForm();

	int def_area(Tuple<int> area_markers);

	void add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym = UNSYM,
	                 int area = ANY, Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ());
        // single equation case
	void add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym = UNSYM,
	                 int area = ANY, Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ())
        {
	  add_matrix_form(0, 0, fn, ord, sym, area, ext);

        }
	void add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, int area = ANY,
	                          Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ());
        // single equation case
	void add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord, int area = ANY,
	                 Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ())
        {
	  add_matrix_form_surf(0, 0, fn, ord, area, ext);

        }

	void add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, int area = ANY, 
                             Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ());
        // single equation case
	void add_vector_form(vector_form_val_t fn, vector_form_ord_t ord, int area = ANY, 
                             Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ())
        {
  	  add_vector_form(0, fn, ord, area, ext);
        };
	void add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, int area = ANY, 
                                  Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ());
        // single equation case
	void add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord, int area = ANY, 
                                  Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ()) 
        {
	  add_vector_form_surf(0, fn, ord, area, ext); 

        };

	void set_ext_fns(void *fn, Tuple<MeshFunction*> ext = Tuple<MeshFunction*> ());

	order3_t get_int_order();
	bool is_matrix_free() { return is_matfree; }

protected:
	int neq;
	bool is_matfree;

	struct Area {
		std::vector<int> markers;
	};

	std::vector<Area> areas;

	struct MatrixFormVol {
		int i, j, sym, area;
		matrix_form_val_t fn; // callback for evaluating the form
		matrix_form_ord_t ord; // callback to determine the integration order
		std::vector<MeshFunction *> ext; // external functions
	};
	struct MatrixFormSurf {
		int i, j, area;
		matrix_form_val_t fn;
		matrix_form_ord_t ord;
		std::vector<MeshFunction *> ext;
	};
	struct VectorFormVol {
		int i, area;
		vector_form_val_t fn;
		vector_form_ord_t ord;
		std::vector<MeshFunction *> ext;
	};
	struct VectorFormSurf {
		int i, area;
		vector_form_val_t fn;
		vector_form_ord_t ord;
		std::vector<MeshFunction *> ext;
	};

	std::vector<MatrixFormVol> mfvol;
	std::vector<MatrixFormSurf> mfsurf;
	std::vector<VectorFormVol> vfvol;
	std::vector<VectorFormSurf> vfsurf;

	struct Stage {
		std::vector<int> idx;
		std::vector<Mesh *> meshes;
		std::vector<Transformable *> fns;
		std::vector<MeshFunction *> ext;

		std::vector<MatrixFormVol *> mfvol;
		std::vector<MatrixFormSurf *> mfsurf;
		std::vector<VectorFormVol *> vfvol;
		std::vector<VectorFormSurf *> vfsurf;

		std::set<int> idx_set;
		std::set<unsigned> seq_set;
		std::set<MeshFunction *> ext_set;
	};

	void get_stages(Space **spaces, std::vector<Stage> &stages, bool rhsonly);
	bool **get_blocks();

	bool is_in_area(int marker, int area) const
	{
		return area >= 0 ? area == marker : is_in_area_2(marker, area);
	}

	bool is_sym() const
	{
		return false; /* not impl. yet */
	}

private:
	Stage *find_stage(std::vector<Stage> &stages, int ii, int jj, Mesh *m1, Mesh *m2,
	                  std::vector<MeshFunction *> &ext);

	bool is_in_area_2(int marker, int area) const;

	// FIXME: pretty dumb to test this in such a way
	bool is_linear() {
		return mfvol.size() > 0 || mfsurf.size() > 0 || vfvol.size() > 0 || vfsurf.size() > 0;
	}

	friend class LinearProblem;
	friend class DiscreteProblem;
	friend class Precond;
};

#endif /* _WEAKFORM_H_ */
