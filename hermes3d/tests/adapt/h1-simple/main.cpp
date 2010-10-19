#define H3D_REPORT_WARN
#define H3D_REPORT_INFO
#define H3D_REPORT_VERBOSE
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

// Test to verify that hp-adaptivity works correctly.
//
// Starts with linear elements and should find the solution (just using p refinements)
//

#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

int P_INIT_X, P_INIT_Y, P_INIT_Z;       // Initial polynomial degree of all mesh elements taken from the command line arguments passed.
const double ERR_STOP = 1.0;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                        // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;           // Adaptivity process stops when the number of degrees of freedom grows
                                        // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
const double TOLERANCE = 0.001;		                // error tolerance in percent
const double THRESHOLD = 0.5;		                  // error threshold for element refinement

int m = 2;
int n = 2;
int o = 2;

// error should be smaller than this epsilon
#define EPS								10e-10F

double fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z) {
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

//

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, Func<res_t> *u_ext[], Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e, ExtData<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, Func<res_t> *u_ext[], Func<f_t> *u, Geom<f_t> *e, ExtData<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) 
{
	int res = ERR_SUCCESS;

	if (argc < 5) error("Not enough parameters.");

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	H3DReader mloader;
	if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);
  
  P_INIT_X = atoi(args[2]);
  P_INIT_Y = atoi(args[3]);
  P_INIT_Z = atoi(args[4]);
	printf("* Setting the space up with initial polynomial degrees %d, %d, %d.\n", P_INIT_X, P_INIT_Y, P_INIT_Z);
	H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM, HERMES_ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

	// Time measurement.
	TimePeriod cpu_time;
	cpu_time.tick();
  
	// Initialize the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
	initialize_solution_environment(matrix_solver, argc, args);

	// Adaptivity loop.
  int as = 1; 
	bool done = false;
	do {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
  	Space* ref_space = construct_refined_space(&space,1 , H3D_H3D_H3D_REFT_HEX_XYZ);
  
    out_orders_vtk(ref_space, "space", as);
	
	  // Initialize discrete problem.
	  bool is_linear = true;
	  DiscreteProblem lp(&wf, ref_space, is_linear);
		
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

    // Assemble the reference problem.
    info("Assembling on reference mesh (ndof: %d).", Space::get_num_dofs(ref_space));
    lp.assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system on reference mesh. If successful, obtain the solution.
    info("Solving on reference mesh.");
    Solution ref_sln(ref_space->get_mesh());
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else {
		  printf("Matrix solver failed.\n");
		  res = ERR_FAILURE;
	  }
    
    // Time measurement.
    cpu_time.tick();

    // Project the reference solution on the coarse mesh.
    Solution sln(space.get_mesh());
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

    // Time measurement.
    cpu_time.tick();

	  // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d.", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    info("err_est_rel: %g%%.", err_est_rel);

	  // If err_est_rel is too large, adapt the mesh. 
    if (err_est_rel < ERR_STOP) {
		  done = true;
      ExactSolution ex_sln(&mesh, exact_solution);
		  
      // Calculation of norms of the solution and the error.
		  double h1_sln_norm = h1_norm(&sln);
		  double h1_err_norm = h1_error(&sln, &ex_sln);
		  printf(" - H1 solution norm: % le\n", h1_sln_norm);
		  printf(" - H1 error norm: % le\n", h1_err_norm);
			
		  double l2_sln_norm = l2_norm(&sln);
		  double l2_err_norm = l2_error(&sln, &ex_sln);
		  printf(" - L2 solution norm: % le\n", l2_sln_norm);
		  printf(" - L2 error norm: % le\n", l2_err_norm);
			
      // If the calculated solution is not precise enough, the returned flag will be failure.
		  if (h1_err_norm > EPS || l2_err_norm > EPS)
		    res = ERR_FAILURE;
	 
      break;
	  }	
    else {
      info("Adapting coarse mesh.");
      adaptivity->adapt(THRESHOLD);
    }

    // If we reached the maximum allowed number of degrees of freedom, set the return flag to failure.
    if (Space::get_num_dofs(&space) >= NDOF_STOP)
    {
      done = true;
      res = ERR_FAILURE;
    }

	  // Clean up.
    delete ref_space->get_mesh();
    delete ref_space;
    delete matrix;
    delete rhs;
    delete solver;
    delete adaptivity;

    // Increase the counter of performed adaptivity steps.
    as++;
	} while (!done);

  // Properly terminate the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
  finalize_solution_environment(matrix_solver);
  
  return res;
}

