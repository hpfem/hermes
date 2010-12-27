//  This benchmark solves a steady state compressible advection-reaction problem with discontinuous solution.
//
//  The following discretization methods may be used:
//    - standard continuous Galerkin method,
//    - streamline upwind Petrov-Galerkin method,
//    - discontinuous Galerkin method with upwind flux,
//    - discontinuous Galerkin method with jump stabilization according to [1]
//  They differ in the mechanism that suppresses unphysical oscillations which develop near the discontinuities.
//
//  PDE: \nabla \cdot (\beta u) + cu = 0, where \beta = ( 10y^2 - 12x + 1, 1 + y ), c = -\nabla\cdot\beta = 11.
//
//  Domain: Square (0,1)x(0,1).
//
//  BC: Dirichlet,  u = 1 on [0,0.5] x {0} and on {0} x [0,0.5]  (bottom left corner),
//                  u = \sin \pi y^2 on {1} x [0,1]  (right side), 
//                  u = 0 elsewhere.
// 
//  Ref.: 
//    [1] F. Brezzi, L. D. Marini, and E. Suli: 
//        Discontinuous Galerkin methods for first-order hyperbolic problems.
//        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.4.333
//    [2] P. Houston, R. Rannacher, E. Süli:
//        A posteriori error analysis for stabilised finite element approximations of transport problems.
//        Comput. Meth. Appl. Mech. Engrg. 190 (2000), pp. 1483-1508.
//    [3] R. Codina:
//        Comparison of some finite element methods for solving the diffusion-convection-reaction equation.
//        Comput. Meth. Appl. Mech. Engrg. 156 (1998), pp. 185-210.
//

#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes2d.h>

// For writing results into file.
#include <fstream>
#include <iterator>

using namespace RefinementSelectors;

enum DiscretizationMethod
{
  CG,
  CG_STAB_SUPG,
  DG_UPWIND,
  DG_JUMP_STAB
};

const std::string method_names[4] = 
{
  "unstabilized continuous Galerkin",
  "streamline upwind Petrov-Galerkin",
  "discontinuous Galerkin - upwind flux",
  "discontinuous Galerkin - Brezzi-Marini-Suli jump stabilization"
};

const DiscretizationMethod method = DG_UPWIND;

const int INIT_REF = 0;                           // Number of initial uniform mesh refinements.
const int P_INIT_H = 0, P_INIT_V = 0;             // Initial polynomial degrees of mesh elements in vertical and horizontal
                                                  // directions.
const double THRESHOLD = 0.20;                    // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = -1 ... do not refine.
                                                  // STRATEGY =  0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY =  1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY =  2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const int ORDER_INCREASE = 1;                                                  
const double ERR_STOP = 1.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

// Diameter of the smallest element of the grid (for SUPG stabilization).
double MIN_DIAM;

// Flow field and reaction term.

template<typename Real>
inline Real fn_a(Real x, Real y) {
  return 10.*sqr(y) - 12.*x + 1.;
}
template<typename Real>
inline Real fn_b(Real x, Real y) {
  return 1. + y;
}
template<typename Real>
inline Real fn_c(Real x, Real y) {
  return 11;  // = -div(beta)
}

// Boundary conditions.

const int BDY_NONZERO_CONSTANT_INFLOW = 1;
const int BDY_ZERO_INFLOW = 2;
const int BDY_VARYING_INFLOW = 3;
const int BDY_OUTFLOW = 4;

template<typename Real, typename Scalar>
inline Scalar fn_g(Real x, Real y)
{
  return sqr(sin(M_PI*y));
}

// Function values for Dirichlet boundary conditions.
template<typename Real, typename Scalar>
Scalar essential_bc_values(int ess_bdy_marker, Real x, Real y)
{
  switch (ess_bdy_marker) 
  {
    case BDY_NONZERO_CONSTANT_INFLOW:
      return 1.0; break;
    case BDY_VARYING_INFLOW:
      return fn_g<Real,Scalar>(x,y); break;
  };
  
  return 0.0;
}

// Forcing source term.
template<typename Real>
Real F(Real x, Real y)
{
  return 0;
}

// Weak forms.
#include "forms.cpp"

// Integral over the specified boundary edge.
double bdry_integral(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  RefMap rm;
    
  Mesh* mesh = sln->get_mesh();
  
  double integral = 0.0;
  Element* el;
  for_all_active_elements(el, mesh)
  {
    for (int e = 0; e < el->nvert; e++)
    {
      if (el->en[e]->bnd && el->en[e]->marker == marker)
      {
        sln->set_active_element(el);
        rm.set_active_element(el);
        
        // Set quadrature order.
        int eo = quad->get_edge_points(e, g_safe_max_order);
        sln->set_quad_order(eo, H2D_FN_VAL);
                
        // Obtain quadrature points, corresponding solution values and tangent vectors (for computing outward normal).
        double* x = rm.get_phys_x(eo);
        double* y = rm.get_phys_y(eo);
        scalar* z = sln->get_fn_values();
        double3* t = rm.get_tangent(e,eo);
        
        // Add contribution from current edge.
        int np = quad->get_num_points(eo);
        double3* pt = quad->get_points(eo);
                
        for (int i = 0; i < np; i++)
        {
          double beta_dot_n = dot2<double>( t[i][1],-t[i][0],fn_a<double>(x[i],y[i]),fn_b<double>(x[i],y[i]) );
          // Weights sum up to two on every edge, therefore the division by two must be present.
          integral += 0.5 * pt[i][2] * t[i][2] * sin(M_PI*x[i]*0.5) * z[i] * beta_dot_n; 
        } 
      }
    }
  }
  
  return integral;
}

int main(int argc, char* args[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Objects storing info about boundary conditions.
  BCTypes bc_types;
  BCValues bc_vals;
  
  // Create a space and refinement selector appropriate for the selected discretization method.
  Space *space;
  ProjBasedSelector *selector;
  ProjNormType norm;
  if (method == CG || method == CG_STAB_SUPG)
  {  
    bc_types.add_bc_dirichlet(Hermes::Tuple<int>(BDY_NONZERO_CONSTANT_INFLOW, BDY_VARYING_INFLOW, BDY_ZERO_INFLOW));
    bc_types.add_bc_neumann(BDY_OUTFLOW);
    bc_vals.add_zero(BDY_ZERO_INFLOW);
    bc_vals.add_const(BDY_NONZERO_CONSTANT_INFLOW, 1.0);
    bc_vals.add_function(BDY_VARYING_INFLOW, fn_g<double>);
    
    space = new H1Space(&mesh, &bc_types, &bc_vals, Ord2(P_INIT_H, P_INIT_V));
    norm = HERMES_L2_NORM;
    //norm = HERMES_H1_NORM;
    
    if (STRATEGY > -1)
    {
      //selector = new H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      selector->set_error_weights(0.5, 2.0, 1.414); // Prefer h-refinement over p-refinement.
    }
  }
  else
  {
    // All boundaries are treated as natural, enforcing the Dirichlet conditions in a weak sense.
    bc_types.add_bc_neumann(Hermes::Tuple<int>(BDY_NONZERO_CONSTANT_INFLOW, BDY_ZERO_INFLOW, BDY_VARYING_INFLOW, BDY_OUTFLOW));
    space = new L2Space(&mesh, &bc_types, Ord2(P_INIT_H, P_INIT_V));
    norm = HERMES_L2_NORM;
    
    if (STRATEGY > -1)
    {
      selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);   
      selector->set_error_weights(0.5, 2.0, 1.414); // Prefer h-refinement over p-refinemet.
    }
  }
  
  // Initialize the weak formulation.
  info("Discretization method: %s", method_names[method].c_str());
    
  WeakForm wf;
  switch(method)
  {
    case CG:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form_surf(callback(cg_biform_surf), BDY_OUTFLOW);
      break;
    case CG_STAB_SUPG:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form_surf(callback(cg_biform_surf), BDY_OUTFLOW);
      wf.add_matrix_form(callback(stabilization_biform_supg));
      break; 
    case DG_UPWIND:
    case DG_JUMP_STAB:
      wf.add_matrix_form(callback(dg_volumetric_biform));
      wf.add_matrix_form_surf(callback(dg_boundary_biform));
      wf.add_vector_form_surf(callback(dg_boundary_liform));
      wf.add_matrix_form_surf(callback(dg_interface_biform), H2D_DG_INNER_EDGE);
      break;
  }      
      
  Solution sln;
  Solution ref_sln;

  // Display the mesh.
  OrderView oview("Distribution of polynomial orders", new WinGeom(0, 0, 500, 400));
  oview.show(space);
  
  ScalarView sview("Solution", new WinGeom(510, 0, 500, 400));
  sview.set_3d_mode(true);
  sview.set_palette(H2DV_PT_HUESCALE);
  sview.fix_scale_width(60);
  //sview.set_min_max_range(0, 1);
  
  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  // Setup data structures for solving the discrete algebraic problem.
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
  
  int as = 1; bool done = false;
  char title[128];
  Space* actual_sln_space;
  do
  {
    info("---- Adaptivity step %d:", as);
    
    if (STRATEGY == -1)
      actual_sln_space = space;
    else
      // Construct globally refined reference mesh and setup reference space.
      actual_sln_space = construct_refined_space(space, ORDER_INCREASE);

    // Find diameter of the smallest element (used by the empirical SUPG stabilization).
    if (method == CG_STAB_SUPG)
    {
      Element *e;
      MIN_DIAM = 1e7;
      for_all_active_elements(e, actual_sln_space->get_mesh())
        if (e->get_diameter() < MIN_DIAM)
          MIN_DIAM = e->get_diameter();
    }
    
    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, actual_sln_space, is_linear);
    
    // Speed-up assembly when 0-th order elements are used.
    bool is_fvm = true;
    Element *e;
    for_all_active_elements(e, actual_sln_space->get_mesh()) 
    {
      int order2d = actual_sln_space->get_element_order(e->id);
      if (H2D_GET_H_ORDER(order2d) > 0 || H2D_GET_V_ORDER(order2d) > 0)
      {
        is_fvm = false;
        break;
      }
    }
    if (is_fvm)
      dp->set_fvm();
    
    // Assemble the linear problem.
    info("Assembling (ndof: %d%s", Space::get_num_dofs(actual_sln_space), is_fvm ? ", FVM)." : ").");
    dp->assemble(matrix, rhs);

    // Solve the linear system. If successful, obtain the solution.
    info("Solving.");
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), actual_sln_space, &ref_sln);
    else error ("Matrix solver failed.\n");
    
    // Skip visualization time.
    cpu_time.tick();
    
    // View the intermediate solution.
    sprintf(title, "Step-%d solution", as);
    sview.set_title(title);
    sview.show(&ref_sln);
    sprintf(title, "Distribution of polynomial orders (%d DOF).", Space::get_num_dofs(actual_sln_space)); 
    oview.set_title(title);
    oview.show(actual_sln_space);
    
    cpu_time.tick(HERMES_SKIP);
    
    // Calculate integral along the outflow boundary (top side of the rectangle).
    double outflow = bdry_integral(&ref_sln, BDY_OUTFLOW);
    info("Integrated outflow: %f, relative error w.r.t. the ref. value (~0.2465): %g%%", 
         outflow, 100*(outflow - 0.246500343856481)/0.246500343856481);
    
    if (STRATEGY == -1) done = true;  // Do not adapt.
    else
    {
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(space, &ref_sln, &sln, matrix_solver, norm);  
          
      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      Adapt* adaptivity = new Adapt(space, norm);
      bool solutions_for_adapt = true;
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                          HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
  
      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
          Space::get_num_dofs(space), Space::get_num_dofs(actual_sln_space), err_est_rel);

      // Skip graphing time.
      cpu_time.tick();

      // Add entry to DOF and CPU convergence graphs.
      graph_dof_est.add_values(Space::get_num_dofs(space), err_est_rel);
      graph_dof_est.save("conv_dof_est.dat");
      graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
      graph_cpu_est.save("conv_cpu_est.dat");
      
      cpu_time.tick(HERMES_SKIP);

      // If err_est_rel too large, adapt the mesh.
      if (err_est_rel < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(space) >= NDOF_STOP) 
        {
          done = true;
          break;
        }
      }

      // Clean up current adaptation step.
      delete adaptivity;
      if(done == false)
      {
        delete actual_sln_space->get_mesh();
        delete actual_sln_space;
      }
      
      as++;
    }
    
    delete dp;
  }
  while (done == false);

  // Destroy the matrix solver.
  delete solver;
  delete matrix;
  delete rhs;
  
  info("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for keyboard or mouse input.
  View::wait();
   
  // Save the final solution in 100x100 points of the domain for comparison with the semi-analytic solution.
  double y = 0.0, step = 0.01;
  int npts = int(1./step+0.5);
  double *res = new double [npts*npts];
  double *p = res;
  std::ofstream fs("out.dat");
  for (int i = 0; i < npts; i++, y+=step)
  {
    double x = 0.0;
    for (int j = 0; j < npts; j++, x+=step)
      *p++ = ref_sln.get_pt_value(x, y);    
  }
  std::copy(res, res+npts*npts, std::ostream_iterator<double>(fs, "\n"));
  
  fs.close();
  delete [] res;
  
  // Destroy spaces and refinement selector object.
  if (STRATEGY > -1)
  {
    delete space;
    delete actual_sln_space->get_mesh();
    delete selector;
  }
  delete actual_sln_space;
  
  return 0;
}
