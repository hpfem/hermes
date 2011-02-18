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

#ifndef _ADAPT_H_
#define _ADAPT_H_

#include "../../../hermes_common/vector.h"
#include "../weakform/weakform.h"

/// hp-Adaptivity module
///
/// TODO
///
///
/// @ingroup hp-adaptivity

// Constants used by Adapt::calc_error().
#define HERMES_TOTAL_ERROR_REL  0x00  
                                ///< A flag which defines interpretation of the total error. \ingroup g_adapt
                                ///  The total error is divided by the norm and therefore it should be in a range [0, 1].
                                ///  \note Used by Adapt::calc_errors_internal().. This flag is mutually exclusive with ::H3D_TOTAL_ERROR_ABS.
#define HERMES_TOTAL_ERROR_ABS  0x01  
                                ///< A flag which defines interpretation of the total error. \ingroup g_adapt
                                ///  The total error is absolute, i.e., it is an integral over squares of differencies.
                                ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::H3D_TOTAL_ERROR_REL.
#define HERMES_ELEMENT_ERROR_REL 0x00 
                                ///< A flag which defines interpretation of an error of an element. \ingroup g_adapt
                                ///  An error of an element is a square of an error divided by a square of a norm of a corresponding component.
                                ///  When norms of 2 components are very different (e.g. microwave heating), it can help.
                                ///  Navier-stokes on different meshes work only when absolute error (see ::H3D_ELEMENT_ERROR_ABS) is used.
                                ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::H3D_ELEMENT_ERROR_ABS.
#define HERMES_ELEMENT_ERROR_ABS 0x10 
                                ///< A flag which defines interpretation of of an error of an element. \ingroup g_adapt
                                ///  An error of an element is a square of an asolute error, i.e., it is an integral over squares of differencies.
                                ///  \note Used by Adapt::calc_errors_internal(). This flag is mutually exclusive with ::H3D_ELEMENT_ERROR_REL.
#define HERMES_TOTAL_ERROR_MASK 0x0F 
                                ///< A mask which mask-out total error type. Used by Adapt::calc_errors_internal(). \internal
#define HERMES_ELEMENT_ERROR_MASK 0xF0 
                                ///< A mask which mask-out element error type. Used by Adapt::calc_errors_internal(). \internal

class HERMES_API Adapt {
public:
	/// Initializes the class.
  /// Constructor. Suitable for problems where various solution components belong to different spaces (L2, H1, Hcurl, 
  /// Hdiv). If proj_norms are not specified, they are expected to be set later by set_error_form.
  Adapt(Hermes::vector<Space *> sp, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>()) {this->init(sp, proj_norms);};
  Adapt(Space* sp, ProjNormType proj_norm = HERMES_H1_NORM) {this->init(Hermes::vector<Space *> (sp), Hermes::vector<ProjNormType> (proj_norm));};
	void init(Hermes::vector<Space *> sp, Hermes::vector<ProjNormType> proj_norms);
	~Adapt();

	typedef
		scalar (*biform_val_t)(int n, double *wt, Func<scalar> *u_ext[], Func<scalar> *u, Func<scalar> *v,
		                       Geom<double> *e, ExtData<scalar> *);
	typedef
		Ord (*biform_ord_t)(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
		                      Geom<Ord> *e, ExtData<Ord> *);

	

	/// Selects elements to refine (based on results from calc_err_est() or calc_energy_error())
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

  /// Type-safe version of calc_err_est() for one solution.
  /// @param[in] solutions_for_adapt - if sln and rsln are the solutions error of which is used in the function adapt().
	double calc_err_est(Solution *sln, Solution *rsln, bool solutions_for_adapt = true, unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS)
	{
		if (num != 1) EXIT("Wrong number of solutions.");
		return calc_err_est(Hermes::vector<Solution *> (sln), Hermes::vector<Solution *> (rsln), solutions_for_adapt, error_flags);
	}

	/// Calculates the error of the solution. 'n' must be the same
	/// as 'num' in the constructor. After that, n coarse solution
	/// pointers are passed, followed by n fine solution pointers.
  /// @param[in] solutions_for_adapt - if slns and rslns are the solutions error of which is used in the function adapt().
	double calc_err_est(Hermes::vector<Solution *> slns, Hermes::vector<Solution *> rslns, bool solutions_for_adapt = true, unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, Hermes::vector<double>* component_errors = NULL)
  {
    return calc_err_internal(slns, rslns, error_flags, component_errors, solutions_for_adapt);
  }

  /// Type-safe version of calc_err_exact() for one solution.
  /// @param[in] solutions_for_adapt - if sln and rsln are the solutions error of which is used in the function adapt().
	double calc_err_exact(Solution *sln, Solution *rsln, bool solutions_for_adapt = true, unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS)
	{
		if (num != 1) EXIT("Wrong number of solutions.");
		return calc_err_exact(Hermes::vector<Solution *> (sln), Hermes::vector<Solution *> (rsln), solutions_for_adapt, error_flags);
	}

	/// Calculates the error of the solution. 'n' must be the same
	/// as 'num' in the constructor. After that, n coarse solution
	/// pointers are passed, followed by n exact solution pointers.
  /// @param[in] solutions_for_adapt - if slns and rslns are the solutions error of which is used in the function adapt().
	double calc_err_exact(Hermes::vector<Solution *> slns, Hermes::vector<Solution *> rslns, bool solutions_for_adapt = true, unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, Hermes::vector<double>* component_errors = NULL)
  {
    return calc_err_internal(slns, rslns, error_flags, component_errors, solutions_for_adapt);
  }

protected:

  /// Internal, calculates the error of the solution. 'n' must be the same
	/// as 'num' in the constructor. After that, first n solution
	/// pointers are passed, followed by second n solution pointers.
	double calc_err_internal(Hermes::vector<Solution *> slns, Hermes::vector<Solution *> rslns, unsigned int error_flags, Hermes::vector<double>* component_errors, bool solutions_for_adapt);

	// parameters
	bool h_only;
	int strategy;
	int max_order;
	bool aniso;
	double exponent;			// exponent used in score formula (see get_optimal_refinement())
  Hermes::vector<ProjNormType> proj_norms;

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
	void get_optimal_refinement(Mesh *mesh, Element *e, const Ord3 &order, Solution *rsln,
	                            Shapeset *ss, int &split, Ord3 p[8]);
	double get_projection_error(Element *e, int split, int son, const Ord3 &order, Solution *rsln, Shapeset *ss);
	int get_dof_count(int split, Ord3 order[]);

	Ord3 get_form_order(int marker, const Ord3 &ordu, const Ord3 &ordv, RefMap *ru,
	                        biform_ord_t bf_ord);
	scalar eval_error(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *sln1,
	                  MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2);
	scalar eval_norm(int marker, biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1,
	                 MeshFunction *rsln2);

	struct ProjKey {
		int split;			// transformation index
		int son;
		Ord3 order;			// element order

    bool operator <(const ProjKey & other) const {
      if(this->split < other.split)
        return true;
      else if(this->split > other.split)
        return false;
      else
        if(this->son < other.son)
          return true;
        else if(this->son > other.son)
          return false;
        else
          if(this->order.order < other.order.order)
            return true;
          else return false;
    };

		ProjKey(int t, int s, const Ord3 &o) {
			split = t;
			son = s;
			order = o;
		}
	};
	std::map<ProjKey, double> proj_err;				// cache for projection errors

	// debugging
	FILE *log_file;
};

#endif
