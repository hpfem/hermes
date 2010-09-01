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

// Solving a linear problem using JFNK method
//
// Dirichlet BC:
//
//   -\Delta u = f in \Omega
//           u = g on d\Omega
//
// Neumann BC:
//
//   -\Delta u + u = f in \Omega
//           du/dn = g on d\Omega
//
// Newton BC:
//
//   -\Delta u = f in \Omega
//   du/dn + u = g on d\Omega
//
//   u(x,y,z) = x^2 + y^2 + z^2
//

#include "config.h"
#include <math.h>
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#include <common/utils.h>

// error should be smaller than this epsilon
#define EPS								10e-10F

// helper macro
#define grad_grad(u, v) (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i])

// functions
template<typename T>
T fnc(T x, T y, T z)
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

#ifdef LIN_DIRICHLET

// case with dirichlet BC

template<typename T>
T dfnc(T x, T y, T z)
{
	return -6.0;
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
res_t form_0(int n, double *wt, fn_t<res_t> *u[], fn_t<f_t> *vi, geom_t<f_t> *e,
             user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi) - dfnc(e->x[i], e->y[i], e->z[i]) * vi->fn[i]);
	return res;
}

// precond
template<typename f_t, typename res_t>
res_t precond_0_0(int n, double *wt, fn_t<f_t> *u[0], fn_t<f_t> *vi, fn_t<f_t> *vj, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#elif defined LIN_NEUMANN

// case with neumann BC

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = 2;
	T ddyy = 2;
	T ddzz = 2;

	return -(ddxx + ddyy + ddzz) + fnc(x, y, z);
}

BCType bc_types(int marker)
{
	return BC_NATURAL;
}

template<typename f_t, typename res_t>
res_t form_0(int n, double *wt, fn_t<f_t> *u[0], fn_t<f_t> *vi, geom_t<f_t> *e,
             user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi)
				+ u[0]->fn[i] * vi->fn[i]
				- dfnc(e->x[i], e->y[i], e->z[i]) * vi->fn[i]);
	return res;
}

template<typename f_t, typename res_t>
res_t form_0_surf(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t result = 0;
	for (int i = 0; i < n; i++) {
		res_t dx = 2 * e->x[i];
		res_t dy = 2 * e->y[i];
		res_t dz = 2 * e->z[i];

		result += wt[i] * (vi->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return -result;
}

// precond
template<typename f_t, typename res_t>
res_t precond_0_0(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, fn_t<f_t> *vj, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#elif defined LIN_NEWTON

// case with newton BC

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = 2;
	T ddyy = 2;
	T ddzz = 2;

	return -(ddxx + ddyy + ddzz);
}

BCType bc_types(int marker)
{
	return BC_NATURAL;
}

template<typename f_t, typename res_t>
res_t form_0(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, geom_t<f_t> *e,
             user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi)
				- dfnc(e->x[i], e->y[i], e->z[i]) * vi->fn[i]);
	return res;
}

template<typename f_t, typename res_t>
res_t form_0_surf(int n, double *wt, fn_t<f_t> *u[], fn_t<f_t> *vi, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0;
	for (int i = 0; i < n; i++) {
		res_t dx = 2 * e->x[i];
		res_t dy = 2 * e->y[i];
		res_t dz = 2 * e->z[i];

		res += wt[i] * (u[0]->fn[i] * vi->fn[i]
				- (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc(e->x[i], e->y[i], e->z[i])) * vi->fn[i]);
	}
	return res;
}

#elif defined NLN_DIRICHLET

// nonlinear case with dirichlet BC

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename T>
inline T f(T x, T y, T z)
{
	return -10 * (x*x + y*y + z*z);
}


template<typename f_t, typename res_t>
res_t form_0(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *vi, geom_t<f_t> *e,
             user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u_ext[0]->fn[i] * grad_grad(u_ext[0], vi) - f(e->x[i], e->y[i], e->z[i]) * vi->fn[i]);
	return res;
}

// precond
template<typename f_t, typename res_t>
res_t precond_0_0(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *vi, fn_t<f_t> *vj, geom_t<f_t> *e,
                  user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#endif


// main ////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int res = ERR_SUCCESS;
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("loading mesh file '%s'\n", argv[1]);

	H1ShapesetLobattoHex shapeset;
	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	space.set_essential_bc_values(essential_bc_values);
#endif

	int mx = 2;
	order3_t order(mx, mx, mx);
	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

	WeakForm wf(true);
	wf.add_resform(form_0<double, scalar>, form_0<ord_t, ord_t>);
#if defined LIN_NEUMANN || defined LIN_NEWTON
	wf.add_resform_surf(form_0_surf<double, scalar>, form_0_surf<ord_t, ord_t>);
#endif
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	// preconditioner
	wf.add_jacform(precond_0_0<double, scalar>, precond_0_0<ord_t, ord_t>, SYM);
#endif

	FeProblem fep(&wf);
	fep.set_space(&space);

#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	// use ML preconditioner to speed-up things
	MlPrecond pc("sa");
	pc.set_param("max levels", 6);
	pc.set_param("increasing or decreasing", "decreasing");
	pc.set_param("aggregation: type", "MIS");
	pc.set_param("coarse: type", "Amesos-KLU");
#endif

	NoxSolver solver(&fep);
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
//	solver.set_precond(&pc);
#endif

	bool solved = solver.solve();
	if (solved) {
		double *s = solver.get_solution();
		Solution sln(&mesh);
		sln.set_fe_solution(&space, s);

		Solution ex_sln(&mesh);
		ex_sln.set_exact(exact_solution);

		double h1_err = h1_error(&sln, &ex_sln);
		printf("  - H1 error norm:      % le\n", h1_err);
		double l2_err = l2_error(&sln, &ex_sln);
		printf("  - L2 error norm:      % le\n", l2_err);

		if (h1_err > EPS || l2_err > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#ifdef OUTPUT_DIR
		// uncomment for debugging purposes
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
			warning("Can not open '%s' for writing.", of_name);
		}
#endif

	}
	else
		res = ERR_FAILURE;

	return res;
}
