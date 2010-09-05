#include "hermes2d.h"
#define H2D_REPORT_FILE "test.log"

/*
 * Benchmark used for testing linear forms in DG.
 */


const int P_INIT = 2;          // Polynomial degree of mesh elements
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

scalar exact_solution(double x, double y, double &dx, double &dy)
{
  dx = 3*sqr(x);
  dy = 3*sqr(y);
  return x*x*x + y*y*y;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += wt[i] * (6*e->x[i]+6*e->y[i]) * v->val[i];
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_boundary(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += -wt[i] * ( 3*sqr(e->x[i])*e->nx[i] + 3*sqr(e->y[i])*e->ny[i] ) * v->val[i];
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_inner(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++){
    result += -wt[i] * 0.5 * ( (ext->fn[0]->get_dx_central(i) + ext->fn[0]->get_dx_neighbor(i))*e->nx[i] + 
                               (ext->fn[0]->get_dy_central(i) + ext->fn[0]->get_dy_neighbor(i))*e->ny[i] ) * v->val[i];
  }
  return result;
}

// boundary conditions
BCType bc_types(int marker)
{
   return BC_NONE;
}

int main(int argc, char* argv[])
{
  // Decide which mesh to refine and how from the command line argument.

  enum ref_type {NO_REF, DEFAULT_REF, TEST_REF};
  ref_type refinement = NO_REF;
  std::string mesh_file = "square.mesh";

  if (argc < 2) {
    info("Using default refinement of square.mesh.");
    refinement = DEFAULT_REF;
  }  else {
    if (!strcmp(argv[1], "test")) {
      info("Using testing refinement of square.mesh.");
      refinement = TEST_REF;
    }
    else {
      mesh_file = argv[1];
      info("Using mesh %s", argv[1]);
    }
  }

  // Load and refine the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(mesh_file.c_str(), &mesh);

  switch (refinement) {
    case DEFAULT_REF:
      // Initial uniform refinement.
      for (int i=0; i<3; i++) mesh.refine_all_elements();

      // Two additional refinements of selected elements.
      mesh.refine_element(29);
      mesh.refine_element(87);

      break;
    case TEST_REF:
      mesh.refine_element(0);
      mesh.refine_element(1);
      mesh.refine_element(7);
      break;
  }
  
  mesh.refine_all_elements();

  // create the L2 space
  L2Space space(&mesh,P_INIT);
  space.set_bc_types(bc_types);

  // display the mesh
  OrderView oview("info_neighbor", 100, 100, 500, 500);
  oview.show(&space);
  BaseView bview("info_neighbor", 100, 100, 500, 500);
  bview.show(&space);
  
  Solution sln;
  Solution xprev(&mesh, &exact_solution);

  // initialize the weak formulation
    WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_vector_form_surf(callback(linear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_vector_form_surf(callback(linear_form_inner), H2D_DG_INNER_EDGE, &xprev);

  // assemble and solve the finite element problem
  solve_linear(&space, &wf, matrix_solver, &sln);

  // visualize the solution
  ScalarView view1("Solution", 600, 100, 500, 500);
  view1.show(&sln);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

