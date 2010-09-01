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

#ifndef _ADAPT_H1_H_
#define _ADAPT_H1_H_

#include "../tuple.h"
#include "../weakform.h"

/// hp-Adaptivity module for H1 space
///
/// TODO
///
///
/// @ingroup hp-adaptivity
class H1Adapt {
public:
	/// Initializes the class. 'num' is the number of mesh-space pairs to be adapted.
	/// After 'num', exactly that many space pointers must follow.
        H1Adapt(Tuple<Space *> sp) {this->init(sp);};
        H1Adapt(Space* sp) {this->init(Tuple<Space *> (sp));};
	void init(Tuple<Space *> sp);
	~H1Adapt();

	typedef
		scalar (*biform_val_t)(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<scalar> *u, fn_t<scalar> *v,
		                       geom_t<double> *e, user_data_t<scalar> *);
	typedef
		ord_t (*biform_ord_t)(int n, double *wt, fn_t<ord_t> *u_ext[], fn_t<ord_t> *u, fn_t<ord_t> *v,
		                      geom_t<ord_t> *e, user_data_t<ord_t> *);

	/// Type-safe version of calc_error_n() for one solution.
	double calc_error(Solution *sln, Solution *rsln)
	{
		if (num != 1) EXIT("Wrong number of solutions.");
		return calc_error_n(Tuple<Solution *> (sln), Tuple<Solution *> (rsln));
	}

	/// Calculates the error of the solution. 'n' must be the same
	/// as 'num' in the constructor. After that, n coarse solution
	/// pointers are passed, followed by n fine solution pointers.
	double calc_error_n(Tuple<Solution *> slns, Tuple<Solution *> rslns);

	/// Selects elements to refine (based on results from calc_error() or calc_energy_error())
	/// and performs their optimal hp-refinement.
	void adapt(double thr);

	/// Get the number of elements refined in the last adaptivity iteration
	int get_num_refined_elements() { return reft_elems; }

	/// Return the time need to calculate the error (in secs)
	double get_error_time() { return error_time; }

	/// Return the time need to the adaptivity step (in secs)
	double get_adapt_time() { return adapt_time; }

	/// Set the type of adaptivity
	/// @param[in] h_only - true if h-adaptivity should be performed
	void set_type(int h_only) { this->h_only = h_only; }

	/// Set the type of adaptivity
	/// @param[in] h_only - true if h-adaptivity should be performed
	void set_exponent(double exp) { this->exponent = exp; }

	/// Allow/disable aniso refinements
	/// @param[in] aniso - set to true to allow aniso refts, false to disable them
	void set_aniso(int aniso) { this->aniso = aniso; }

	/// Set strategy type
	/// @param[in] strat - type of the strategy
	/// - 0 = adapt until the error drop below certain threshold
	/// - 1 = adapt until specified percentage of elements is refined
	void set_strategy(int strat) { strategy = strat; }

	/// Set the log file
	/// - use for debugging the adaptivity
	void set_log_file(FILE *log_file) { this->log_file = log_file; }

	/// Set the bilinear form on position [i,j]
	/// TODO: be more verbose here
	void set_error_form(int i, int j, biform_val_t biform, biform_ord_t biform_ord);

	/// Set maximal order to be used
	void set_max_order(int order)
	{
		if (order <= H3D_MAX_ELEMENT_ORDER) max_order = order;
		else error("Cannot set order greater than max elem. order (%d)", H3D_MAX_ELEMENT_ORDER);
	}

protected:
	// parameters
	bool h_only;
	int strategy;
	int max_order;
	bool aniso;
	double exponent;			// exponent used in score formula (see get_optimal_refinement())

	// spaces & solutions
	int num;					// number of spaces to work with
	Space **spaces;
	Solution **sln;
	Solution **rsln;

	// element error arrays
	double **errors;
	double  *norms;
	bool    have_errors;
	double  total_err;
	int2 *esort;
	int   nact;
	int reft_elems;

	double error_time;			// time needed to calculate error
	double adapt_time;			// time needed for adaptivity step

	// bilinear forms to calculate error
	biform_val_t **form;
	biform_ord_t **ord;

	/// Used by adapt(). Can be utilized in specialized adaptivity
	/// procedures, for which adapt() is not sufficient.
	void get_optimal_refinement(Mesh *mesh, Element *e, const order3_t &order, Solution *rsln,
	                            Shapeset *ss, int &split, order3_t p[8]);
	double get_projection_error(Element *e, int split, int son, const order3_t &order, Solution *rsln, Shapeset *ss);
	int get_dof_count(int split, order3_t order[]);

	order3_t get_form_order(int marker, const order3_t &ordu, const order3_t &ordv, RefMap *ru,
	                        biform_ord_t bf_ord);
	scalar eval_error(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *sln1,
	                  MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2);
	scalar eval_norm(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1,
	                 MeshFunction *rsln2);

	struct ProjKey {
		int split;			// transformation index
		int son;
		order3_t order;			// element order

		ProjKey(int t, int s, const order3_t &o) {
			split = t;
			son = s;
			order = o;
		}
	};
	Map<ProjKey, double> proj_err;				// cache for projection errors

	// debugging
	FILE *log_file;
};

#endif
