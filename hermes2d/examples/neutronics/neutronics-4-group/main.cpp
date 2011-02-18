#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. The corresponding diffusion operator is given by (r = x, z = y):
//
//	\nabla \cdot D \nabla \phi = \frac{1}{x} (x D \phi_x)_x  + (D \phi_y)_y 
//
// BC:
//
// Homogeneous neumann on symmetry axis,
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 3;                       // Number of initial uniform mesh refinements.
const int P_INIT_1 = 4,                           // Initial polynomial degree for approximation of group 1 fluxes.
          P_INIT_2 = 4,                           // Initial polynomial degree for approximation of group 2 fluxes.
          P_INIT_3 = 4,                           // Initial polynomial degree for approximation of group 3 fluxes.
          P_INIT_4 = 4;                           // Initial polynomial degree for approximation of group 4 fluxes.
const double ERROR_STOP = 1e-5;                   // Tolerance for the eigenvalue.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)

// Initial eigenvalue approximation.
double k_eff = 1.0;         

// Element markers.
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary markers.
const int BDY_VACUUM = 1;
const int BDY_SYM = 2;

// Physical data of the problem for the 4 energy groups.
#include "physical_parameters.cpp"
// Weak forms.
#include "forms.cpp"

// Source function.
void source_fn(int n, Hermes::vector<scalar*> values, scalar* out)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = (nu[1][0] * Sf[1][0] * values.at(0)[i] +
    nu[1][1] * Sf[1][1] * values.at(1)[i] +
    nu[1][2] * Sf[1][2] * values.at(2)[i] +
    nu[1][3] * Sf[1][3] * values.at(3)[i]);
  }
}

// Integral over the active core.
double integrate(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }

  return 2.0 * M_PI * integral;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Solution variables.
  Solution sln1, sln2, sln3, sln4;
  Solution iter1, iter2, iter3, iter4;
  Hermes::vector<Solution*> solutions(&sln1, &sln2, &sln3, &sln4);

  // Define initial conditions.
  info("Setting initial conditions.");
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(BDY_SYM);
  bc_types.add_bc_newton(BDY_VACUUM);

  // Create H1 spaces with default shapesets.
  H1Space space1(&mesh, &bc_types, P_INIT_1);
  H1Space space2(&mesh, &bc_types, P_INIT_2);
  H1Space space3(&mesh, &bc_types, P_INIT_3);
  H1Space space4(&mesh, &bc_types, P_INIT_4);
  Hermes::vector<Space*> spaces(&space1, &space2, &space3, &space4);
  
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&space1, &space2, &space3, &space4));
  info("ndof = %d.", ndof);
  
  // Initialize views.
  ScalarView view1("Neutron flux 1", new WinGeom(0, 0, 320, 600));
  ScalarView view2("Neutron flux 2", new WinGeom(350, 0, 320, 600));
  ScalarView view3("Neutron flux 3", new WinGeom(700, 0, 320, 600));
  ScalarView view4("Neutron flux 4", new WinGeom(1050, 0, 320, 600));
  
  // Do not show meshes.
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);
  
  // Initialize the weak formulation.
  WeakForm wf(4);
  wf.add_matrix_form(0, 0, callback(biform_0_0), HERMES_SYM);
  wf.add_matrix_form(1, 1, callback(biform_1_1), HERMES_SYM);
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(2, 2, callback(biform_2_2), HERMES_SYM);
  wf.add_matrix_form(2, 1, callback(biform_2_1));
  wf.add_matrix_form(3, 3, callback(biform_3_3), HERMES_SYM);
  wf.add_matrix_form(3, 2, callback(biform_3_2));
  wf.add_vector_form(0, callback(liform_0), marker_core, Hermes::vector<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(1, callback(liform_1), marker_core, Hermes::vector<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(2, callback(liform_2), marker_core, Hermes::vector<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(3, callback(liform_3), marker_core, Hermes::vector<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), BDY_VACUUM);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), BDY_VACUUM);
  wf.add_matrix_form_surf(2, 2, callback(biform_surf_2_2), BDY_VACUUM);
  wf.add_matrix_form_surf(3, 3, callback(biform_surf_3_3), BDY_VACUUM);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, spaces, is_linear);
  
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Time measurement.
  TimePeriod cpu_time, solver_time;
  
  // Main power iteration loop:
  int iter = 1; bool done = false;
  bool rhs_only = false;
  
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
  do
  {
    info("------------ Power iteration %d:", iter);
    
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs, rhs_only);
    
    /* 
    // Testing the factorization reuse schemes for direct solvers.
    if (iter == 10)  
      solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING);
    if (iter == 20)
      solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING_AND_SCALING);
    if (iter == 30) 
      solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
    */
 
    info("Solving the matrix problem by %s.", MatrixSolverNames[matrix_solver].c_str());
    solver_time.tick(HERMES_SKIP);  
    bool solved = solver->solve();  
    solver_time.tick();
    
    if(solved)
      Solution::vector_to_solutions(solver->get_solution(), spaces, solutions);
    else
      error ("Matrix solver failed.\n");

    // Show intermediate solutions.
    // view1.show(&sln1);    
    // view2.show(&sln2);
    // view3.show(&sln3);    
    // view4.show(&sln4);

    SimpleFilter source(source_fn, Hermes::vector<MeshFunction*>(&sln1, &sln2, &sln3, &sln4));
    SimpleFilter source_prev(source_fn, Hermes::vector<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));

    // Compute eigenvalue.
    double k_new = k_eff * (integrate(&source, marker_core) / integrate(&source_prev, marker_core));
    info("Largest eigenvalue: %.8g, rel. difference from previous it.: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;

    // Update eigenvalue.
    k_eff = k_new;
    
    if (!done)
    {
      // Save solutions for the next iteration.
      iter1.copy(&sln1);    
      iter2.copy(&sln2);
      iter3.copy(&sln3);    
      iter4.copy(&sln4);
      
      // Don't need to reassemble the system matrix in further iterations,
      // only the rhs changes to reflect the progressively updated source.
      rhs_only = true;

      iter++;
    }
  }
  while (!done);
  
  // Time measurement.
  cpu_time.tick();
  solver_time.tick(HERMES_SKIP);
  
  // Print timing information.
  verbose("Average solver time for one power iteration: %g s", solver_time.accumulated() / iter);
  
  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  // Show solutions.
  view1.show(&sln1);
  view2.show(&sln2);
  view3.show(&sln3);    
  view4.show(&sln4);
  
  // Skip visualization time.
  cpu_time.tick(HERMES_SKIP);

  // Print timing information.
  verbose("Total running time: %g s", cpu_time.accumulated());
    
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
