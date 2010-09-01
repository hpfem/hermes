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

#ifndef _NOX_SOLVER_H_
#define _NOX_SOLVER_H_

#include "../solver.h"
#include "../discrete_problem.h"
#include "../precond.h"
#include "epetra.h"

#ifdef HAVE_NOX
#include <NOX.H>
#include <NOX_Epetra.H>
#endif

class NoxProblemInterface;

/// Encapsulation of NOX nonlinear solver
///
/// @ingroup solvers
class NoxSolver : public Solver {
public:
	NoxSolver(DiscreteProblem *problem);
	virtual ~NoxSolver();

	bool set_init_sln(double *ic);
	bool set_init_sln(EpetraVector *ic);
	virtual bool solve();

	int get_num_iters() { return num_iters; }

	// settings for the solver
	void set_nl_method(const char *par);
	void set_output_info(int flags) { output_flags = flags; }

	// linear solver setters
	void set_ls_type(const char *type) { ls_type = type; }
	void set_ls_max_iters(int iters) { ls_max_iters = iters; }
	void set_ls_tolerance(double tolerance) { ls_tolerance = tolerance; }
	void set_ls_sizeof_krylov_subspace(int size) { ls_sizeof_krylov_subspace = size; }

	// convergence params
	void set_conv_iters(int iters)        { conv.max_iters = iters; }
	void set_conv_abs_resid(double resid) { conv_flag.absresid = 1; conv.abs_resid = resid; }
	void set_conv_rel_resid(double resid) { conv_flag.relresid = 1; conv.rel_resid = resid; }
	void set_conv_update(double update)   { conv_flag.update = 1; conv.update = update; }
	void set_conv_wrms(double rtol, double atol) {
		conv_flag.wrms = 1;
		conv.wrms_rtol = rtol;
		conv.wrms_atol = atol;
	}

	void set_precond(Precond *pc);

	double get_assembly_time();
	double get_precond_time();

protected:
#ifdef HAVE_NOX
	Teuchos::RCP<NoxProblemInterface> interface;
#endif
	int num_iters;
	const char *nl_dir;
	int output_flags;

	const char *ls_type;
	int ls_max_iters;
	double ls_tolerance;
	int ls_sizeof_krylov_subspace;
	// convergence params
	struct conv_t {
		int max_iters;
		double abs_resid;
		double rel_resid;
		double update;
		double wrms_rtol;
		double wrms_atol;
	} conv;

	struct conv_flag_t {
		unsigned absresid:1;
		unsigned relresid:1;
		unsigned wrms:1;
		unsigned update:1;
	} conv_flag;
};

#endif
