#define HERMES_REPORT_ALL
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

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT_1 = 1,                           // Initial polynomial degree for approximation of group 1 fluxes.
          P_INIT_2 = 1,                           // Initial polynomial degree for approximation of group 2 fluxes.
          P_INIT_3 = 2,                           // Initial polynomial degree for approximation of group 3 fluxes.
          P_INIT_4 = 2;                           // Initial polynomial degree for approximation of group 4 fluxes.
const double ERROR_STOP = 1e-5;                   // Tolerance for the eigenvalue.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)

// Initial eigenvalue approximation.
double k_eff = 1.0;         

// Element markers.
std::string reflector = "reflector";
std::string core = "core";

// Boundary markers.
std::string bdy_vacuum = "vacuum boundary";
std::string bdy_symmetry = "symmetry plane";

// Weak forms, input data and some other utility functions.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Solution variables.
  Solution sln1, sln2, sln3, sln4;
  Hermes::vector<Solution*> solutions(&sln1, &sln2, &sln3, &sln4);
  
  // Define initial conditions.
  info("Setting initial conditions.");
  Solution iter1, iter2, iter3, iter4;
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);
  Hermes::vector<MeshFunction*> iterates(&iter1, &iter2, &iter3, &iter4);

  // Create H1 spaces with default shapesets.
  H1Space space1(&mesh, P_INIT_1);
  H1Space space2(&mesh, P_INIT_2);
  H1Space space3(&mesh, P_INIT_3);
  H1Space space4(&mesh, P_INIT_4);
  Hermes::vector<Space*> spaces(&space1, &space2, &space3, &space4);
  
  int ndof = Space::get_num_dofs(spaces);
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
  
  // Load physical data of the problem for the 4 energy groups.
  MaterialPropertyMaps matprop(4);
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_scattering_multigroup_structure(Ss_nnz);
  matprop.set_fission_multigroup_structure(chi_nnz);
  matprop.set_Sigma_a(Sa);
  matprop.set_Sigma_f(Sf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.validate();
  
  std::cout << matprop;
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, iterates, k_eff, bdy_vacuum);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, spaces);
  
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
  
  // Initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  
  // Force the Jacobian assembling in the first iteration.
  bool Jacobian_changed = true;
  
  // In the following iterations, Jacobian will not be changing; its LU factorization
  // may be reused.
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
  
  // Main power iteration loop:
  int it = 1; bool done = false;
  do
  {
    info("------------ Power iteration %d:", it);
    
    info("Newton's method (matrix problem solved by %s).", MatrixSolverNames[matrix_solver].c_str());
    
    memset(coeff_vec, 0.0, ndof*sizeof(scalar)); //TODO: Why it doesn't work without zeroing coeff_vec in each iteration?
    
    solver_time.tick(HERMES_SKIP);      
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, Jacobian_changed, 1e-8, 10, true)) 
      error("Newton's iteration failed.");
    solver_time.tick();
    
    Solution::vector_to_solutions(solver->get_solution(), spaces, solutions);
    
    // Show intermediate solutions.
    view1.show(&sln1);    
    view2.show(&sln2);
    view3.show(&sln3);    
    view4.show(&sln4);
    
    // Compute eigenvalue.
    
    SourceFilter source(solutions, matprop);
    SourceFilter source_prev(iterates, matprop);
    
    double k_new = k_eff * (integrate(&source, core) / integrate(&source_prev, core));
    info("Largest eigenvalue: %.8g, rel. difference from previous it.: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;

    // Update eigenvalue.
    k_eff = k_new;
    wf.update_keff(k_eff);
    
    if (!done)
    {
      // Save solutions for the next iteration.
      iter1.copy(&sln1);    
      iter2.copy(&sln2);
      iter3.copy(&sln3);    
      iter4.copy(&sln4);
      
      // Don't need to reassemble the system matrix in further iterations,
      // only the rhs changes to reflect the progressively updated source.
      Jacobian_changed = false;

      it++;
    }
  }
  while (!done);
  
  delete [] coeff_vec;
  
  // Time measurement.
  cpu_time.tick();
  solver_time.tick(HERMES_SKIP);
  
  // Print timing information.
  verbose("Average solver time for one power iteration: %g s", solver_time.accumulated() / it);
  
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
