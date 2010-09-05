#include "hermes2d.h"
#define H2D_REPORT_FILE "test.log"

/*
 *  Linear transport equation.
 */

// General settings.

const int INIT_REF = 1;       // Number of initial uniform mesh refinements.
const int P_INIT = 3;         // Polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Flux stuff.

template<typename Real>
inline Real calculate_a_dot_v(Real x, Real y, Real vx, Real vy) 
{
  Real norm = std::max(1e-12, sqrt(sqr(x) + sqr(y)));
  return -y/norm*vx + x/norm*vy;
}

template<>
inline Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) 
{
  return Ord(10);
}

template<typename Real, typename Scalar>
inline Scalar upwind_flux(Real u_cent, Real u_neib, Real a_dot_n)
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib); 
}

template<>
inline Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n)
{
  return a_dot_n * (u_cent + u_neib); 
}


// Boundary conditions.

BCType bc_types(int marker)
{
  return BC_NONE;
}
// function values for Dirichlet boundary conditions.
template<typename Real, typename Scalar>
Scalar g(int ess_bdy_marker, Real x, Real y)
{
  if (ess_bdy_marker == 1) return 1; else return 0;
}

template<typename Real>
Real F(Real x, Real y)
{
  return 0;
}

// Weak forms.

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += -wt[i] * u->val[i] * calculate_a_dot_v<Real>(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real,Scalar>(n, wt, F<Real>, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_interface(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real a_dot_n = calculate_a_dot_v<Real>(e->x[i], e->y[i], e->nx[i], e->ny[i]);
    Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
    result += wt[i] * upwind_flux<Real,Scalar>(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_boundary(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = calculate_a_dot_v<Real>(x, y, e->nx[i], e->ny[i]);
    result += wt[i] * upwind_flux<Real,Scalar>(u->val[i], 0, a_dot_n) * v->val[i];
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_boundary(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = calculate_a_dot_v<Real>(x, y, e->nx[i], e->ny[i]);
    result += -wt[i] * upwind_flux<Real,Scalar>(0, g<Real,Scalar>(e->marker,x,y), a_dot_n) * v->val[i];
  }
  
  return result;
}

// MAIN

int main(int argc, char* argv[])
{
  // Load and refine the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial uniform refinement.
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
  
  // create the L2 space
  L2Space space(&mesh,P_INIT);
  space.set_bc_types(bc_types);

  // display the mesh
  OrderView oview("Distribution of polynomial orders", 100, 100, 500, 500);
  oview.show(&space);
  BaseView bview("Distribution of polynomial orders", 600, 100, 500, 500);
  bview.show(&space);
  MeshView mview("Mesh", 100, 600, 500, 500);
  mview.show(&mesh);

  Solution sln;

  // initialize the weak formulation
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_vector_form_surf(callback(linear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_matrix_form_surf(callback(bilinear_form_interface), H2D_DG_INNER_EDGE);

  // assemble and solve the finite element problem
  solve_linear(&space, &wf, matrix_solver, &sln);
  
  // visualize the solution
  ScalarView view1("Solution", 600, 600, 500, 500);
  view1.show(&sln);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

