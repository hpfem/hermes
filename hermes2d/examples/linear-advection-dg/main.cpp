#define H2D_REPORT_INFO
#include "hermes2d.h"

//  This example solves a linear advection equation using Discontinuous Galerkin (DG) method.
//  Its purpose is to show how to evaluate surface matrix forms that take basis functions defined
//  on different elements.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC:		Dirichlet,  u = g where \Beta(x) \cdot n(x) < 0; g = 1 on [0,0.5] x {0}, g = 0 anywhere else.
//				
//  The following parameters can be changed:

const int P_INIT = 0;                             // Polynomial degree of mesh elements.
const int INIT_REF = 1;                           // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Flux definition.
#include "fluxes.cpp"


// Boundary conditions.
BCType bc_types(int marker)
{
  return BC_NONE;
}
// Essential (Dirichlet) boundary condition values.
template<typename Real, typename Scalar>
Scalar g(int ess_bdy_marker, Real x, Real y)
{
  if (ess_bdy_marker == 1) return 1; else return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinement.
  mesh.refine_element(0, 1);
  mesh.refine_element(2);
  mesh.refine_element(3);
  mesh.refine_element(5, 2);
  mesh.refine_element(6); 
  mesh.refine_element(7, 1);
  mesh.refine_element(8, 1); 
  mesh.refine_element(14);
  mesh.refine_element(15);
  mesh.refine_element(16);
  mesh.refine_element(17);
  mesh.refine_element(19);
  mesh.refine_element(20);
  mesh.refine_element(42);
  mesh.refine_element(43);
  mesh.refine_element(44);
  mesh.refine_element(45);
  mesh.refine_element(13, 1);
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  //mesh.refine_all_elements();
  //mesh.refine_all_elements();
  
  // Create an L2 space with default shapeset.
  L2Space space(&mesh,P_INIT);
  space.set_bc_types(bc_types);
  info("ndof: %d", get_num_dofs(&space));

  // Display the mesh.
  OrderView oview("Mesh", new WinGeom(0, 0, 440, 350));
  oview.show(&space);
  BaseView bview("Basis Functions", new WinGeom(450, 0, 400, 350));
  bview.show(&space);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_vector_form_surf(callback(linear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_matrix_form_surf(callback(bilinear_form_interface), H2D_DG_INNER_EDGE);

  // Assemble and solve the finite element problem.
  Solution sln;
  info("Solving...");
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(&wf, &space, is_linear);
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  dp->assemble(matrix, rhs);
  
  // Time measurement.
  //cpu_time.tick();
  /*
  MumpsMatrix* mm = static_cast<MumpsMatrix*>(matrix);
  MumpsVector* mv = static_cast<MumpsVector*>(rhs);
  FILE *file;
  file = fopen("test.txt", "wt");
  mm->dump(file, "A", DF_PLAIN_ASCII);
  mv->dump(file, "x", DF_PLAIN_ASCII);
  fclose(file);
  */
  // Solve the linear system of the reference problem. 
  // If successful, obtain the solution.
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");
  
  // Visualize the solution.
  ScalarView view1("Solution", new WinGeom(860, 0, 400, 350));
  view1.show(&sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

