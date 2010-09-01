// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

//  Testing nonlinear solver using Newton's method

#include "config.h"
#include <math.h>
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#include <common/utils.h>

// error should be smaller than this epsilon
#define EPS								10e-10F

#define grad_grad(u, v) (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i])


#if defined LINEAR

//  Solving a linear problem using nonlinear solver with Newton's method
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity
//  -\Laplace u = f
//
//  Exact solution: u(x,y,z) = x^2 + y^2 + z^2

double fnc(double x, double y, double z)
{
	return x*x + y*y + z*z;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t jacobi_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *vi, fn_t<f_t> *vj, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, vi, vj, e);
}

template<typename f_t, typename res_t>
res_t resid_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *vi, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u_ext[0], vi) + 6.0 * vi->fn[i]);
	return res;
}

#elif defined NONLIN1

//  Solving a nonlinear problem using Newton's method
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity
//  - div[lambda(u) grad u] = f
//                lambda(u) = 1 + u^2
//
//  BC:  T = 100 on the left, top and bottom edges
//       dT/dn = 0 on the face
//
//  Exact solution: u(x,y,z) = 100

// thermal conductivity (temperature-dependent)
// for any u, this function has to be  positive in the entire domain!
template<typename T>
inline T lambda(T temp)  { return 10 + 0.1 * pow(temp, 2); }
// derivate of lambda wrt temperature
template<typename T>
inline T dlambda(T temp) { return 0.2 * temp; }

const int marker_right = 2;

BCType bc_types(int marker)
{
	if (marker == marker_right) return BC_NATURAL;
	else return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return 100.0;
}

template<typename f_t, typename res_t>
res_t jacobi_form(int n, double *wt, fn_t<res_t> *itr[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] *
			(dlambda(itr[0]->fn[i]) * u->fn[i] * grad_grad(itr[0], v) +
			  lambda(itr[0]->fn[i]) * grad_grad(u, v));
	return res;
}

template<typename f_t, typename res_t>
res_t resid_form(int n, double *wt, fn_t<res_t> *itr[], fn_t<f_t> *v, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (lambda(itr[0]->fn[i]) * grad_grad(itr[0], v));
	return res;
}

#elif defined NONLIN2

//  Solving a nonlinear problem using Newton's method
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity
//  - div [u grad u] = f
//
//  BC:  u = x^2 + y^2 + z^2 on dOmega
//
//  Exact solution: u(x,y,z) = x^2 + y^2 + z^2

double fnc(double x, double y, double z)
{
	return x*x + y*y + z*z;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

// BC

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

// forms

template<typename T>
inline T f(T x, T y, T z)
{
	return -10.0 * (x*x + y*y + z*z);
}


template<typename f_t, typename res_t>
res_t jacobi_form(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, fn_t<f_t> *vj, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (vj->fn[i] * grad_grad(u[0], vi) + u[0]->fn[i] * grad_grad(vj, vi));
	return res;
}


template<typename f_t, typename res_t>
res_t resid_form(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u[0]->fn[i] * grad_grad(u[0], vi) - f(e->x[i], e->y[i], e->z[i]) * vi->fn[i]);
	return res;
}

// PROJ

template<typename f_t, typename res_t>
res_t biproj_form(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u->fn[i] * v->fn[i] + grad_grad(u, v));
	return res;
}

template<typename f_t, typename res_t>
res_t liproj_form(int n, double *wt, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	res_t res = 0.0;
	return res;
}

#endif

// main ////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	int res = ERR_SUCCESS;
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("loading mesh file '%s'\n", argv[1]);

	H1ShapesetLobattoHex shapeset;

#if defined NONLIN1
	order3_t order(1, 1, 1);
#else
	order3_t order(2, 2, 2);
#endif
	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

#if defined NONLIN2
	// do L2 projection of zero function
	WeakForm proj_wf;
	proj_wf.add_biform(biproj_form<double, scalar>, biproj_form<ord_t, ord_t>, SYM);
	proj_wf.add_liform(liproj_form<double, scalar>, liproj_form<ord_t, ord_t>);

	LinProblem lp(&proj_wf);
	lp.set_space(&space);

#ifdef WITH_UMFPACK
	UMFPackMatrix m;
	UMFPackVector v;
	UMFPackLinearSolver sl(&m, &v);
#elif defined WITH_MUMPS
	MumpsMatrix m;
	MumpsVector v;
	MumpsSolver sl(&m, &v);
#endif
	lp.assemble(&m, &v);
	sl.solve();

	double *ps = sl.get_solution();
#endif

	printf("* Calculating a solution\n");

	WeakForm wf(1);
	wf.add_jacform(0, 0, jacobi_form<double, scalar>, jacobi_form<ord_t, ord_t>, UNSYM);
	wf.add_resform(0, resid_form<double, scalar>, resid_form<ord_t, ord_t>);

	FeProblem fep(&wf);
	fep.set_spaces(1, &space);

	NoxSolver solver(&fep);
#if defined NONLIN2
	solver.set_init_sln(ps);
#endif
	solver.set_conv_iters(10);

	printf("  - solving..."); fflush(stdout);
	Timer solve_timer;
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	if (solved) {
		printf(" done in %s (%lf secs), iters = %d\n", solve_timer.get_human_time(),
		       solve_timer.get_seconds(), solver.get_num_iters());

		double *s = solver.get_solution();
		Solution sln(&mesh);
		sln.set_fe_solution(&space, s);

		Solution ex_sln(&mesh);
#ifdef NONLIN1
		ex_sln.set_const(100.0);
#else
		ex_sln.set_exact(exact_solution);
#endif
		double h1_err = h1_error(&sln, &ex_sln);
		printf("  - H1 error norm:      % le\n", h1_err);
		double l2_err = l2_error(&sln, &ex_sln);
		printf("  - L2 error norm:      % le\n", l2_err);

		if (h1_err > EPS || l2_err > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}
#ifdef OUTPUT_DIR
		printf("* Output\n");
		// output
		const char *of_name = OUTPUT_DIR "/solution.vtk";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			VtkOutputEngine output(ofile);
			output.out(&sln, "Uh", FN_VAL_0);
			fclose(ofile);
		}
		else {
			warning("Cann not open '%s' for writing.", of_name);
		}

#endif
	}
	else
		res = ERR_FAILURE;

	return res;
}
