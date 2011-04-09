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

#include "../function.h"
#include "forms.h"
#include "../../../hermes_common/vector.h"

// Bilinear form symmetry flag, see WeakForm::add_matrix_form
enum SymFlag {
        HERMES_ANTISYM = -1,
	HERMES_NONSYM = 0,
	HERMES_SYM = 1
};

/// Matrix and vector forms.
typedef scalar (*matrix_form_val_t)(int n, double *wt, Func<scalar> *u_ext[], Func<double> *vi,
	        Func<double> *vj, Geom<double> *e, ExtData<scalar> *);
typedef Ord (*matrix_form_ord_t)(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *vi,
	       Func<Ord> *vj, Geom<Ord> *e, ExtData<Ord> *);
typedef scalar (*vector_form_val_t)(int n, double *wt, Func<scalar> *u_ext[], Func<double> *vi,
	        Geom<double> *e, ExtData<scalar> *);
typedef Ord (*vector_form_ord_t)(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *vi,
	       Geom<Ord> *e, ExtData<Ord> *);

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
class HERMES_API WeakForm {

public:
	WeakForm(int neq = 1, bool mat_free = false);
	virtual ~WeakForm();

	int def_area(Hermes::vector<int> area_markers);

	void add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym = HERMES_NONSYM,
	                 int area = HERMES_ANY_INT, Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ());
        // single equation case
	void add_matrix_form(matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym = HERMES_NONSYM,
	                 int area = HERMES_ANY_INT, Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ())
        {
	  add_matrix_form(0, 0, fn, ord, sym, area, ext);

        }
	void add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, int area = HERMES_ANY_INT,
	                          Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ());
        // single equation case
	void add_matrix_form_surf(matrix_form_val_t fn, matrix_form_ord_t ord, int area = HERMES_ANY_INT,
	                 Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ())
        {
	  add_matrix_form_surf(0, 0, fn, ord, area, ext);

        }

	void add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, int area = HERMES_ANY_INT, 
                             Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ());
        // single equation case
	void add_vector_form(vector_form_val_t fn, vector_form_ord_t ord, int area = HERMES_ANY_INT, 
                             Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ())
        {
  	  add_vector_form(0, fn, ord, area, ext);
        };
	void add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, int area = HERMES_ANY_INT, 
                                  Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ());
        // single equation case
	void add_vector_form_surf(vector_form_val_t fn, vector_form_ord_t ord, int area = HERMES_ANY_INT, 
                                  Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ()) 
        {
	  add_vector_form_surf(0, fn, ord, area, ext); 

        };

	void set_ext_fns(void *fn, Hermes::vector<MeshFunction*> ext = Hermes::vector<MeshFunction*> ());

        /// Returns the number of equations
        int get_neq() { return neq; }

        /// Internal. Used by DiscreteProblem to detect changes in the weakform.
        int get_seq() const { return seq; }

	Ord3 get_int_order();
	bool is_matrix_free() { return is_matfree; }

        int neq;

protected:

        int seq;
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


	// OLD CODE: 
        // void get_stages(Space **spaces, std::vector<Stage> &stages, bool rhsonly);
        void get_stages(Hermes::vector< Space* > spaces, Hermes::vector< Solution* >& u_ext, 
                        std::vector< WeakForm::Stage >& stages, bool want_matrix, bool want_vector);
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

        Stage* find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                          Mesh* m1, Mesh* m2, 
                          std::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext);

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
