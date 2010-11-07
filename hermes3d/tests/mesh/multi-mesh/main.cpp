#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <math.h>
#include <hermes3d.h>

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
// multi-mesh.cc
//
// -Laplace u1 + u2      = f1
//      -Laplace u2 + u3 = f2
//           -Laplace u3 = f3
//
// u1 = x^2 + y^2 + z^2
// u2 = (1-x^2)(1 - y^2)(1- z^2)
// u3 = x

// error should be smaller than this epsilon
#define EPS								10e-10F

#ifdef RHS2

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

// Boundary condition types.
BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                    ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (data->ext[0].fn[i] * v->fn[i]);
	return res;
}

#elif defined SYS

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
//	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return x;
//	return 3.0;
//	return (1 - x*x) * x*x * (1 - y*y) * y*y * (1 - z*z) * z*z;
}

// needed for calculation norms and used by visualizator
double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = -2 * x * (1 - y*y) * (1 - z*z);
//	dy = -2 * (1 - x*x) * y * (1 - z*z);
//	dz = -2 * (1 - x*x) * (1 - y*y) * z;
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = 2 * x * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	dy = 2 * x*x * (1 - x*x) * y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y*y * z*z * (1 - z*z);
//	dz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z*z;
	dx = 1;
	dy = 0;
	dz = 0;

	return u2(x, y, z);
}

// Boundary condition types.
BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

// Boundary condition types.
BCType bc_types_2(int marker)
{
	if (marker == 3) return BC_NATURAL;
	else return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z) {
	return u2(x, y, z);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                        ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                        ExtData<Scalar> *data)
{
//	return 0.0;
	return int_u_v<Real, Scalar>(n, wt, u, v, e);
}

//template<typename T>
//T f1(T x, T y, T z)
//{
//	T ddxx = -2 * (1 - y*y) * (1 - z*z);
//	T ddyy = -2 * (1 - x*x) * (1 - z*z);
//	T ddzz = -2 * (1 - x*x) * (1 - y*y);
//
//	return -(ddxx + ddyy + ddzz) + u2(x, y, z);
//}

template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data)
{
//	return int_F_v<Real, Scalar>(n, wt, f1, u, e);
//	return -3.0 * int_u<Real, Scalar>(n, wt, u, e);
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * ((-6.0 + (e->x[i]) * u->fn[i]);
	return res;
//	int_u<Real, Scalar>(n, wt, u, e);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//template<typename Real, typename Scalar>
//Scalar bilinear_form_2_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
//                        ExtData<Scalar> *data)
//{
//	return int_u_v<Real, Scalar>(n, wt, u, v, e);
//}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                        ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

//template<typename T>
//T f2(T x, T y, T z)
//{
//	T ddxx = 2 * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	T ddyy = 2 * x*x * (1 - x*x) * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * z*z * (1 - z*z);
//	T ddzz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z;
//
//	return -(ddxx + ddyy + ddzz) + u1(x, y, z);
//}

template<typename Real, typename Scalar>
Scalar linear_form_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data)
{
//	return int_F_v<Real, Scalar>(n, wt, f2, u, e);
//	return -6.0 * int_u<Real, Scalar>(n, wt, u, e);
	return 0;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_surf(int np, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data)
{
	Scalar result = 0;
	for (int i = 0; i < np; i++) {
		Scalar dx = 1;
		Scalar dy = 0;
		Scalar dz = 0;

		result += wt[i] * (u->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return result;
}

#elif defined SYS3

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u3(T x, T y, T z)
{
	return x;
}

double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = -2 * x * (1 - y*y) * (1 - z*z);
	dy = -2 * (1 - x*x) * y * (1 - z*z);
	dz = -2 * (1 - x*x) * (1 - y*y) * z;

	return u2(x, y, z);
}

double exact_sln_fn_3(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 1;
	dy = 0;
	dz = 0;

	return u3(x, y, z);
}

// Boundary condition types.
BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

// Boundary condition types.
BCType bc_types_2(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z)
{
	return 0;
}

// Boundary condition types.
BCType bc_types_3(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values_3(int ess_bdy_marker, double x, double y, double z)
{
	return x;
}

// 1. eqn ------------------------------------------------------------------------------------------

template<typename Real, typename Scalar>
Scalar biform_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_1_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	return int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename T>
T rhs1(T x, T y, T z)
{
	return -6.0 + u2(x, y, z);
}

template<typename Real, typename Scalar>
Scalar liform_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (rhs1(e->x[i], e->y[i], e->z[i]) * v->fn[i]);
	return res;
}

// 2. eqn ------------------------------------------------------------------------------------------

template<typename Real, typename Scalar>
Scalar biform_2_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_2_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	return int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename T>
T rhs2(T x, T y, T z)
{
	T laplace = 2 * (1 - y*y) * (1 - z*z) + 2 * (1 - x*x) * (1 - z*z) + 2 * (1 - x*x) * (1 - y*y);
	return laplace + u3(x, y, z);
}

template<typename Real, typename Scalar>
Scalar liform_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data)
{
	return int_F_v<Real, Scalar>(n, wt, rhs2, v, e);
}

// 3. eqn ------------------------------------------------------------------------------------------

template<typename Real, typename Scalar>
Scalar biform_3_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                        ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

// Weak forms. 
#include "forms.cpp"

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 2) error("Not enough parameters.");

	info("* Loading mesh '%s'.", args[1]);
	Mesh mesh1;
	H3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh1)) error("Loading mesh file '%s'.", args[1]);

#if defined RHS2

	Ord3 order(2, 2, 2);
	info("  - Setting uniform order to (%d, %d, %d).", order.x, order.y, order.z);

	info("* Setting the space up.");

	H1Space space(&mesh1, bc_types, essential_bc_values, order);

	int ndofs = space.assign_dofs();
	info("  - Number of DOFs: %d.", ndofs);

	info("* Calculating a solution.");

	// duplicate the mesh
	Mesh mesh2;
	mesh2.copy(mesh1);
	// do some changes
	mesh2.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
	mesh2.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

	Solution fsln(&mesh2);
	fsln.set_const(-6.0);
#else
	// duplicate the mesh
	Mesh mesh2;
	mesh2.copy(mesh1);

	Mesh mesh3;
	mesh3.copy(mesh1);

	// change meshes
	mesh1.refine_all_elements(H3D_REFT_HEX_X);
	mesh2.refine_all_elements(H3D_REFT_HEX_Y);
	mesh3.refine_all_elements(H3D_REFT_HEX_Z);

	info("* Setup spaces.");
	H1Space space1(&mesh1, &shapeset);
	space1.set_bc_types(bc_types_1);
	space1.set_essential_bc_values(essential_bc_values_1);

	Ord3 o1(2, 2, 2);
	info("  - Setting uniform order to (%d, %d, %d).", o1.x, o1.y, o1.z);
	space1.set_uniform_order(o1);

	H1Space space2(&mesh2, &shapeset);
	space2.set_bc_types(bc_types_2);
	space2.set_essential_bc_values(essential_bc_values_2);

	Ord3 o2(2, 2, 2);
	info("  - Setting uniform order to (%d, %d, %d).", o2.x, o2.y, o2.z);
	space2.set_uniform_order(o2);

	H1Space space3(&mesh3, &shapeset);
	space3.set_bc_types(bc_types_3);
	space3.set_essential_bc_values(essential_bc_values_3);

	Ord3 o3(1, 1, 1);
	info("  - Setting uniform order to (%d, %d, %d).", o3.x, o3.y, o3.z);
	space3.set_uniform_order(o3);

	int ndofs = 0;
	ndofs += space1.assign_dofs();
	ndofs += space2.assign_dofs(ndofs);
	ndofs += space3.assign_dofs(ndofs);
	info("  - Number of DOFs: %d.", ndofs);
#endif

#if defined WITH_UMFPACK
	MatrixSolverType matrix_solver = SOLVER_UMFPACK; 
#elif defined WITH_PARDISO
	MatrixSolverType matrix_solver = SOLVER_PARDISO; 
#elif defined WITH_PETSC
	MatrixSolverType matrix_solver = SOLVER_PETSC; 
#elif defined WITH_MUMPS
	MatrixSolverType matrix_solver = SOLVER_MUMPS; 
#endif

#ifdef RHS2
	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &fsln);

	// Initialize discrete problem.
	bool is_linear = true;
	DiscreteProblem dp(&wf, &space, is_linear);
#elif defined SYS3
	WeakForm wf(3);
	wf.add_matrix_form(0, 0, biform_1_1<double, scalar>, biform_1_1<Ord, Ord>, HERMES_SYM);
	wf.add_matrix_form(0, 1, biform_1_2<double, scalar>, biform_1_2<Ord, Ord>, HERMES_UNSYM);
	wf.add_vector_form(0, liform_1<double, scalar>, liform_1<Ord, Ord>);

	wf.add_matrix_form(1, 1, biform_2_2<double, scalar>, biform_2_2<Ord, Ord>, HERMES_SYM);
	wf.add_matrix_form(1, 2, biform_2_3<double, scalar>, biform_2_3<Ord, Ord>, HERMES_UNSYM);
	wf.add_vector_form(1, liform_2<double, scalar>, liform_2<Ord, Ord>);

	wf.add_matrix_form(2, 2, biform_3_3<double, scalar>, biform_3_3<Ord, Ord>, HERMES_SYM);

	// Initialize discrete problem.
	bool is_linear = true;
	DiscreteProblem dp(&wf, Tuple<Space *>(&space1, &space2, &space3), is_linear);
#endif
	// Time measurement.
	TimePeriod cpu_time;
	cpu_time.tick();

	// Initialize the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
	initialize_solution_environment(matrix_solver, argc, args);
  
	// Set up the solver, matrix, and rhs according to the solver selection.
	SparseMatrix* matrix = create_matrix(matrix_solver);
	Vector* rhs = create_vector(matrix_solver);
	Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

	// Initialize the preconditioner in the case of SOLVER_AZTECOO.
	if (matrix_solver == SOLVER_AZTECOO) 
	{
		((AztecOOSolver*) solver)->set_solver(iterative_method);
		((AztecOOSolver*) solver)->set_precond(preconditioner);
		// Using default iteration parameters (see solver/aztecoo.h).
	}

	// Assemble stiffness matrix and load vector.
	dp.assemble(matrix, rhs);

	// Solve the linear system. If successful, obtain the solution.
	info("Solving the linear problem.");
	bool solved = solver->solve();

	// Time measurement.
	cpu_time.tick();
	// Print timing information.
	info("Solution and mesh with polynomial orders saved. Total running time: %g s", cpu_time.accumulated());

	// assemble stiffness matrix
	info("  - assembling... "); fflush(stdout);
	Timer assemble_timer;
	assemble_timer.start();
	dp.assemble(&mat, &rhs);
	assemble_timer.stop();
	info("%s (%lf secs).", assemble_timer.get_human_time(), assemble_timer.get_seconds());

	// solve the stiffness matrix
	info("  - solving... "); fflush(stdout);
	Timer solve_timer;
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();
	info("%s (%lf secs).", solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
#ifdef RHS2
		Timer sln_pre_tmr;
		Solution sln(&mesh1);
		sln_pre_tmr.start();
		sln.set_coeff_vector(&space, solver.get_solution());
		sln_pre_tmr.stop();

		info("* Solution:.");

		// Set exact solution.
		ExactSolution ex_sln(&mesh1, exact_solution);

		// Norm.
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);
		info("  - H1 solution norm:   % le.", h1_sln_norm);
		info("  - H1 error norm:      % le.", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &ex_sln);
		info("  - L2 solution norm:   % le.", l2_sln_norm);
		info("  - L2 error norm:      % le.", l2_err_norm);

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}
#elif defined SYS3
		// Solution 1.
		Solution sln1(&mesh1);
		Solution sln2(&mesh2);
		Solution sln3(&mesh3);

		Solution::vector_to_solution(solver->get_solution(), &space1, &sln1);
		Solution::vector_to_solution(solver->get_solution(), &space2, &sln2);
		Solution::vector_to_solution(solver->get_solution(), &space3, &sln3);

		ExactSolution esln1(&mesh1, exact_sln_fn_1);
		ExactSolution esln2(&mesh2, exact_sln_fn_2);
		ExactSolution esln3(&mesh3, exact_sln_fn_3);

		// Norm.
		double h1_err_norm1 = h1_error(&sln1, &esln1);
		double h1_err_norm2 = h1_error(&sln2, &esln2);
		double h1_err_norm3 = h1_error(&sln3, &esln3);

		double l2_err_norm1 = l2_error(&sln1, &esln1);
		double l2_err_norm2 = l2_error(&sln2, &esln2);
		double l2_err_norm3 = l2_error(&sln3, &esln3);

		info("  - H1 error norm:      % le.", h1_err_norm1);
		info("  - L2 error norm:      % le.", l2_err_norm1);
		if (h1_err_norm1 > EPS || l2_err_norm1 > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}

		info("  - H1 error norm:      % le.", h1_err_norm2);
		info("  - L2 error norm:      % le.", l2_err_norm2);
		if (h1_err_norm2 > EPS || l2_err_norm2 > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}

		info("  - H1 error norm:      % le.", h1_err_norm3);
		info("  - L2 error norm:      % le.", l2_err_norm3);
		if (h1_err_norm3 > EPS || l2_err_norm3 > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}
#endif

#ifdef RHS2
		out_fn_vtk(&sln, "solution");
#elif defined SYS3
		out_fn_vtk(&sln1, "sln1");
		out_fn_vtk(&sln2, "sln2");
		out_fn_vtk(&sln3, "sln3");
#endif
	}
	else
		res = ERR_FAILURE;

	// Print timing information.
	info("Solution and mesh with polynomial orders saved. Total running time: %g s", sln_time.accumulated());

	// Clean up.
	delete matrix;
	delete rhs;
	delete solver;

	// Properly terminate the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
	finalize_solution_environment(matrix_solver);

#ifdef WITH_PETSC
	mat.free();
	rhs.free();

	PetscFinalize();
#endif

	return res;
}
