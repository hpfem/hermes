#include "config.h"
#include <hermes3d.h>
#include "../../../../hermes_common/trace.h"
#include "../../../../hermes_common/error.h"
#include <float.h>

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
// Testing projections
//

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).


#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

#define EPS											10e-10

//#define X2_Y2_Z2_025
//#define X2_Y2_Z2
//#define X3_Y3_Z3
#define XN_YM_ZO

int m = 2, n = 2, o = 2;

double fnc(double x, double y, double z) {
#if defined X2_Y2_Z2_025
  return pow(x*x + y*y + z*z, .25);
#elif defined X2_Y2_Z2
  return x*x + y*y + z*z;
#elif defined X3_Y3_Z3
  return x*x*x + y*y*y + z*z*z;
#elif defined XN_YM_ZO
  return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#endif
}

template<typename T>
T dfnc(T x, T y, T z) {
#if defined X2_Y2_Z2_025
  if ((x*x + y*y + z*z) == 0.0) return -DBL_MAX;
  else
    return -0.75 * pow(x*x + y*y + z*z, -0.75);
#elif defined X3_Y3_Z3
    return -6.0 * x - 6.0 * y - 6.0 * z;
#elif defined X2_Y2_Z2
    return -6.0;
#elif defined XN_YM_ZO
    T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
    T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
    T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
    return -(ddxx + ddyy + ddzz);
#endif
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#if defined X2_Y2_Z2_025
  if ((x*x + y*y + z*z) != 0.0) {
    dx = 0.5 * x * pow(x*x + y*y + z*z, -.75);
    dy = 0.5 * y * pow(x*x + y*y + z*z, -.75);
    dz = 0.5 * z * pow(x*x + y*y + z*z, -.75);
  }
  else {
    // pow(x*x + y*y + z*z, -.75) is not defined
    dx = 0.0;
    dy = 0.0;
    dz = 0.0;
  }
#elif defined X2_Y2_Z2
  dx = 2 * x;
  dy = 2 * y;
  dz = 2 * z;
#elif defined X3_Y3_Z3
  dx = 3 * x*x;
  dy = 3 * y*y;
  dz = 3 * z*z;
#elif defined XN_YM_ZO
  dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
  dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
  dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#endif

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


//

#define CHECK_ERROR \
  printf("Elem #%u: error = % lf\n", e->id, error); \
  if (error > EPS) { \
    ret = ERR_FAILURE; \
    break; \
  }

//
// main
//

int main(int argc, char *args[])
{
  _F_
  int ret = ERROR_SUCCESS;

  if (argc < 3) {
    fprintf(stderr, "ERROR: not enough parameters\n");
    return ERR_FAILURE;
  }

  if (strcmp(args[1], "h1") != 0 && strcmp(args[1], "h1-ipol")) {
    fprintf(stderr, "ERROR: unknown type of the projection\n");
    return ERR_FAILURE;
  }

  Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[2], &mesh)) {
    fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[2]);
    return ERR_FAILURE;
  }

  unsigned int ne = mesh.elements.count();
  // make the mesh for the ref. solution
  mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

#if defined X2_Y2_Z2
  Ord3 o(2, 2, 2);
#elif defined X3_Y3_Z3
  Ord3 o(3, 3, 3);
#elif defined XN_YM_ZO
  Ord3 o(2, 3, 4);
#endif

  H1Space space(&mesh, bc_types, essential_bc_values, o);

  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM, HERMES_ANY);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

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

  // assemble the stiffness matrix
  dp.assemble(matrix, rhs);

  // solve the stiffness matrix
  solver->solve();

  Solution sln(&mesh);
  Solution::vector_to_solution(solver->get_solution(), &space, &sln);

  for (unsigned int idx = mesh.elements.first(); idx <= ne; idx = mesh.elements.next(idx)) {
    Element *e = mesh.elements[idx];

    Ord3 order(4, 4, 4);
    double error;

    Projection *proj;
    if (strcmp(args[1], "h1") == 0) proj = new H1Projection(&sln, e, space.get_shapeset());
    else if (strcmp(args[1], "h1-ipol") == 0) proj = new H1ProjectionIpol(&sln, e, space.get_shapeset());
    else return ERR_FAILURE;

    //
    error = 0.0;
    error += proj->get_error(H3D_REFT_HEX_NONE, -1, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_REFT_HEX_X, 20, order);
    error += proj->get_error(H3D_REFT_HEX_X, 21, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_REFT_HEX_Y, 22, order);
    error += proj->get_error(H3D_REFT_HEX_Y, 23, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_REFT_HEX_Z, 24, order);
    error += proj->get_error(H3D_REFT_HEX_Z, 25, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_H3D_REFT_HEX_XY,  8, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XY,  9, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XY, 10, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XY, 11, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 12, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 13, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 14, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 15, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 16, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 17, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 18, order);
    error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 19, order);
    error = sqrt(error);
    CHECK_ERROR;

    //
    error = 0.0;
    for (int j = 0; j < 8; j++)
      error += proj->get_error(H3D_H3D_H3D_REFT_HEX_XYZ, j, order);
    error = sqrt(error);
    CHECK_ERROR;

    delete proj;
  }

  // Properly terminate the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
  finalize_solution_environment(matrix_solver);

  return ret;
}
