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
//
// error should be smaller than this epsilon

#define EPS								10e-10F

const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;                           // Initial polynomial degree of all mesh elements.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).



// Weak forms. 
#include "definitions.cpp"

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	if (argc < 2) error("Not enough parameters.");

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh1;
	H3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh1)) error("Loading mesh file '%s'\n", args[1]);

#if defined RHS2

	Ord3 order(P_INIT_X, P_INIT_Y, P_INIT_Z);
	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	
	// Create an H1 space with default shapeset.
	printf("* Setting the space up\n");
	H1Space space(&mesh1, bc_types, essential_bc_values, order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

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

	printf("* Setup spaces\n");
	Ord3 o1(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o1.x, o1.y, o1.z);
	H1Space space1(&mesh1, bc_types_1, essential_bc_values_1, o1);

	Ord3 o2(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o2.x, o2.y, o2.z);
	H1Space space2(&mesh2, bc_types_2, essential_bc_values_2, o2);

	Ord3 o3(1, 1, 1);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o3.x, o3.y, o3.z);
	H1Space space3(&mesh3, bc_types_3, essential_bc_values_3, o3);

	int ndofs = 0;
	ndofs += space1.assign_dofs();
	ndofs += space2.assign_dofs(ndofs);
	ndofs += space3.assign_dofs(ndofs);
	printf("  - Number of DOFs: %d\n", ndofs);
#endif

#if defined WITH_UMFPACK
	MatrixSolverType matrix_solver = SOLVER_UMFPACK; 
#elif defined WITH_PETSC
	MatrixSolverType matrix_solver = SOLVER_PETSC; 
#elif defined WITH_MUMPS
	MatrixSolverType matrix_solver = SOLVER_MUMPS; 
#endif

#ifdef RHS2
	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY_INT, &fsln);

	// Initialize discrete problem.
	bool is_linear = true;
	DiscreteProblem dp(&wf, &space, is_linear);
#elif defined SYS3
	WeakForm wf(3);
	wf.add_matrix_form(0, 0, biform_1_1<double, scalar>, biform_1_1<Ord, Ord>, HERMES_SYM);
	wf.add_matrix_form(0, 1, biform_1_2<double, scalar>, biform_1_2<Ord, Ord>, HERMES_NONSYM);
	wf.add_vector_form(0, liform_1<double, scalar>, liform_1<Ord, Ord>);

	wf.add_matrix_form(1, 1, biform_2_2<double, scalar>, biform_2_2<Ord, Ord>, HERMES_SYM);
	wf.add_matrix_form(1, 2, biform_2_3<double, scalar>, biform_2_3<Ord, Ord>, HERMES_NONSYM);
	wf.add_vector_form(1, liform_2<double, scalar>, liform_2<Ord, Ord>);

	wf.add_matrix_form(2, 2, biform_3_3<double, scalar>, biform_3_3<Ord, Ord>, HERMES_SYM);

	// Initialize discrete problem.
	bool is_linear = true;
	DiscreteProblem dp(&wf, Hermes::vector<Space *>(&space1, &space2, &space3), is_linear);
#endif
	// Time measurement.
	TimePeriod cpu_time;
	cpu_time.tick();
  
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

	// Time measurement.
	TimePeriod sln_time;
	sln_time.tick();

	if (solved) {
#ifdef RHS2
		// Solve the linear system. If successful, obtain the solution.
		info("Solving the linear problem.");
                Solution sln(&mesh1);
		Solution::vector_to_solution(solver->get_solution(), &space, &sln);

		// Set exact solution.
		ExactSolution ex_sln(&mesh1, exact_solution);

		// Norm.
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);
		printf("  - H1 solution norm:   % le\n", h1_sln_norm);
		printf("  - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &ex_sln);
		printf("  - L2 solution norm:   % le\n", l2_sln_norm);
		printf("  - L2 error norm:      % le\n", l2_err_norm);

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

		printf("  - H1 error norm:      % le\n", h1_err_norm1);
		printf("  - L2 error norm:      % le\n", l2_err_norm1);
		if (h1_err_norm1 > EPS || l2_err_norm1 > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}

		printf("  - H1 error norm:      % le\n", h1_err_norm2);
		printf("  - L2 error norm:      % le\n", l2_err_norm2);
		if (h1_err_norm2 > EPS || l2_err_norm2 > EPS) {
			// Calculated solution is not enough precise.
			res = ERR_FAILURE;
		}

		printf("  - H1 error norm:      % le\n", h1_err_norm3);
		printf("  - L2 error norm:      % le\n", l2_err_norm3);
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

	return res;
}
