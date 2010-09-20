#define H2D_REPORT_INFO
#include "hermes2d.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//	It is intended to show how evalutation of surface matrix forms that take basis functions defined
//	on different elements work.
//  It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC:		Dirichlet,  u = g where \Beta(x) \cdot n(x) < 0, where g = 1 on [0,0.5] x {0}, g = anywhere else.
//				
//  The following parameters can be changed:

using namespace RefinementSelectors;

const int INIT_REF = 2;                         // Number of initial uniform mesh refinements.
const int P_INIT = 0;                           // Polynomial degree of mesh elements.

const double THRESHOLD = 0.2;                   // This is a quantitative parameter of the adapt(...) function and
                                                // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                         // Adaptive strategy:
                                                // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                //   error is processed. If more elements have similar errors, refine
                                                //   all to keep the mesh symmetric.
                                                // STRATEGY = 1 ... refine all elements whose error is larger
                                                //   than THRESHOLD times maximum element error.
                                                // STRATEGY = 2 ... refine all elements whose error is larger
                                                //   than THRESHOLD.
                                                // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;         // Predefined list of element refinement candidates. Possible values are
                                                // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                 // Maximum allowed level of hanging nodes:
                                                // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                // Note that regular meshes are not supported, this is due to
                                                // their notoriously bad performance.
const double ERR_STOP = 3.0;                    // Stopping criterion for adaptivity (rel. error tolerance between the
                                                // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                    // Default value is 1.0. This parameter influences the selection of
                                                // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                    // Adaptivity process stops when the number of degrees of freedom grows
                                                // over this limit. This is to prevent h-adaptivity to go on forever.


MatrixSolverType matrix_solver = SOLVER_UMFPACK;// Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                // SOLVER_MUMPS, and more are coming.

// Flux definition.

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
// Function values for Dirichlet boundary conditions.
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
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Create the L2 space.
  L2Space space(&mesh,P_INIT);
  space.set_bc_types(bc_types);

  // Initialize refinement selector.
  L2ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  
  // Disable weighting of refinement candidates.
  selector.set_error_weights(1, 1, 1);

  SimpleGraph graph_dof_est, graph_cpu_est;
  TimePeriod cpu_time;
  cpu_time.tick();

  // Display the mesh.
  OrderView oview("Coarse mesh", 100, 100, 500, 500);
  oview.show(&space);
  BaseView bview("Distribution of polynomial orders", 600, 100, 500, 500);
  bview.show(&space);

  Solution sln;
	Solution ref_sln;

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_vector_form_surf(callback(linear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_matrix_form_surf(callback(bilinear_form_interface), H2D_DG_INNER_EDGE);

  ScalarView view1("Solution", 600, 600, 500, 500);

  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space.get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space.dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(&space, order_increase);

    // Solve the reference problem.
    solve_linear(ref_space, &wf, matrix_solver, &ref_sln);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(&space, H2D_L2_NORM, &ref_sln, &sln, NULL); 

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    view1.show(&sln);
    oview.show(&space);
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(&space, H2D_L2_NORM);
    hp.set_solutions(&sln, &ref_sln);
    double err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
 
    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", get_num_dofs(&space), get_num_dofs(ref_space), err_est_rel);

       // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(&space) >= NDOF_STOP) done = true;
    }

    // Free reference space and mesh.
    ref_space->free();
		delete ref_mesh;

    as++;
  }
  while (done == false);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

