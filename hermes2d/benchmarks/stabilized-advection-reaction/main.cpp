//  This benchmark solves a steady state compressible advection-reaction problem with discontinuous solution.
//
//  The following discretization methods may be used:
//    - standard continuous Galerkin method,
//    - streamline upwind Petrov-Galerkin method,
//    - discontinuous Galerkin method with upwind flux.
//  They differ in the mechanism that suppresses unphysical oscillations which develop near the discontinuities.
//
//  PDE: \beta\cdot\nabla u + cu = f, where \beta = ( 10y^2 - 12x + 1, 1 + y ), c = f = 0.
//
//  Domain: Square (0,1)x(0,1).
//
//  BC: Dirichlet,  u = 1 on [0,0.5] x {0} and on {0} x [0,0.5]  (bottom left corner),
//                  u = \sin \pi y^2 on {1} x [0,1]  (right side), 
//                  u = 0 elsewhere.
// 
//  Author: Milan Hanus (University of West Bohemia, Pilsen, Czech Republic).
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

#include <iostream>

// For writing results into file.
#include <fstream>
#include <iterator>

using namespace RefinementSelectors;

enum DiscretizationMethod
{
  CG,
  SUPG,
  DG
};

const std::string method_names[4] = 
{
  "unstabilized continuous Galerkin",
  "streamline upwind Petrov-Galerkin",
  "discontinuous Galerkin",
};

const DiscretizationMethod method = DG;
const bool SAVE_FINAL_SOLUTION = false;           // Save the final solution at specified points for comparison with the
                                                  // semi-analytic solution in Mathematica?
                                                  
const int INIT_REF = 0;                           // Number of initial uniform mesh refinements.
const int P_INIT = 0;                             // Initial polynomial degrees of mesh elements.
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
const CandList CAND_LIST = H2D_HP_ISO;            // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const int ORDER_INCREASE = 1;                     // Difference in polynomial orders of the coarse and the reference spaces.
const double ERR_STOP = 0.0001;                   // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 8000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

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
  if (method == DG)
    return 11.; // -div(beta)
  else
    return 0.; 
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
  return 0.0;
}

// Weak forms.
#include "forms.cpp"

// Exact solution.
#include "exact.h"

// Weighting function for the target functional.
double weight_fn(double x, double y)
{
  return sin(M_PI*x*0.5);
}

// Weighted integral over the specified boundary edge.
double bdry_integral(MeshFunction* sln, int marker, double (*weight_fn_ptr)(double,double))
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
          integral += 0.5 * pt[i][2] * t[i][2] * (*weight_fn_ptr)(x[i],y[i]) * z[i] * beta_dot_n; 
        } 
      }
    }
  }
  
  return integral;
}

// Construct a string representing the program options.
void make_str_from_program_options(std::stringstream& str)
{ 
  switch (method)
  {
    case CG:
      str << "cg";
      break;
    case SUPG:
      str << "supg";
      break;
    case DG:
      str << "dg";
      break;
  }
  
  if (STRATEGY > -1)
  {
    switch (CAND_LIST) 
    {
      case H2D_H_ANISO:
      case H2D_H_ISO:
        str << "_h" << P_INIT;
        break;
      case H2D_P_ANISO:
      case H2D_P_ISO:
        str << "_p" << INIT_REF;
        break;
      default:
        str << "_hp";
        break;
    }
    switch (CAND_LIST) 
    {
      case H2D_H_ANISO:
      case H2D_P_ANISO:
      case H2D_HP_ANISO:
        str << "_aniso";
        break;
      case H2D_H_ISO:
      case H2D_P_ISO:
      case H2D_HP_ISO:
        str << "_iso";
        break;
      case H2D_HP_ANISO_H:
        str << "_anisoh";
        break;
      case H2D_HP_ANISO_P:
        str << "_anisop";
        break;
    }
  }
  else
  {
    str << "_h-" << INIT_REF << "_p-" << P_INIT;
  }
}

int main(int argc, char* args[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Object representing types of condtions prescribed on each part of the domain boundary.
  BCTypes bc_types;
  
  // Create a space and refinement selector appropriate for the selected discretization method.
  Space *space;
  ProjBasedSelector *selector;
  ProjNormType norm;
  
  // For both methods, all inflow boundaries are treated as natural, enforcing the Dirichlet conditions in a weak sense. 
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_NONZERO_CONSTANT_INFLOW, BDY_ZERO_INFLOW, BDY_VARYING_INFLOW));
  if (method != DG)
  { 
    space = new H1Space(&mesh, &bc_types, Ord2(P_INIT, P_INIT));
    norm = HERMES_L2_NORM;
    //norm = HERMES_H1_NORM;
    
    if (STRATEGY > -1)
    {
      selector = new H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      selector->set_error_weights(1.15, 1.0, 1.414); // Prefer h-refinement over p-refinement.
      //selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      //selector->set_error_weights(0.5, 2.0, 1.414); // Prefer h-refinement over p-refinement.
    }
  }
  else
  {
    // For DGM, there is also a bilinear form on the outflow boundary.
    bc_types.add_bc_neumann(BDY_OUTFLOW);
    space = new L2Space(&mesh, &bc_types, Ord2(P_INIT, P_INIT));
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
  wf.add_vector_form(callback(source_liform));
  switch(method)
  {
    case CG:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form_surf(callback(cg_boundary_biform));
      wf.add_vector_form_surf(callback(cg_boundary_liform));
      break;
    case SUPG:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form_surf(callback(cg_boundary_biform));
      wf.add_vector_form_surf(callback(cg_boundary_liform));
      wf.add_matrix_form(callback(stabilization_biform_supg));
      break; 
    case DG:
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
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_ex, graph_cpu_ex, graph_dof_outfl, graph_cpu_outfl;
  
  // Time measurement.
  TimePeriod cpu_time, clk_time;
  cpu_time.tick();
  clk_time.tick();
  
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
  
  // Load the exact solution evaluated at the Gauss-Kronrod quadrature points.
  SemiAnalyticSolution exact_sln("exact/sol_GaussKronrod50.map");
  
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
    
    // Process the intermediate solution, but don't accumulate cpu time.
    cpu_time.tick();
    
    // View the intermediate solution.
    sprintf(title, "Step-%d solution", as);
    sview.set_title(title);
    sview.show(&ref_sln);
    sprintf(title, "Distribution of polynomial orders (%d DOF).", Space::get_num_dofs(actual_sln_space)); 
    oview.set_title(title);
    oview.show(actual_sln_space);
    
    // Calculate relative error w.r.t. exact solution.
    info("Calculating relative L2 error w.r.t. exact solution.");
    double err_exact_rel = exact_sln.get_l2_rel_err(&ref_sln);
    info("Relative L2 error w.r.t. exact solution: %g%%", err_exact_rel*100);
    
    // Calculate integral along the outflow boundary (top side of the rectangle).
    double outflow = bdry_integral(&ref_sln, BDY_OUTFLOW, &weight_fn);
    double err_outflow = std::abs(outflow - 0.246500343856481)/0.246500343856481;
    info("Integrated outflow: %f, relative error w.r.t. the ref. value (~0.2465): %g%%", 
         outflow, 100*err_outflow);
    
    cpu_time.tick(HERMES_SKIP);
         
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
                          HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);
                          
      int ndof_fine = Space::get_num_dofs(actual_sln_space);
      int ndof_coarse = Space::get_num_dofs(space);
  
      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
           ndof_coarse, ndof_fine, err_est_rel*100);

      // Get current time, use it in convergence graphs, and skip the graphing time.
      cpu_time.tick();
      double t = cpu_time.accumulated();

      // Add entry to DOF and CPU convergence graphs.
      graph_dof_est.add_values(ndof_fine, err_est_rel);
      graph_cpu_est.add_values(t, err_est_rel);
      graph_dof_ex.add_values(ndof_fine, err_exact_rel);
      graph_cpu_ex.add_values(t, err_exact_rel);
      graph_dof_outfl.add_values(ndof_fine, err_outflow);
      graph_cpu_outfl.add_values(t, err_outflow);
      
      cpu_time.tick(HERMES_SKIP);

      // If err_est_rel too large, adapt the mesh.
      if (err_est_rel*100 < ERR_STOP) done = true;
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
        as++;
      }
    }
    
    delete dp;
  }
  while (done == false);

  // Destroy the matrix solver.
  delete solver;
  delete matrix;
  delete rhs;
  
  clk_time.tick();
  info("Total running time: %g s", clk_time.accumulated());
  info("Total cpu time: %g s", cpu_time.accumulated());
  
  // Wait for keyboard or mouse input.
  View::wait();
  
  // Get the string summarizing the currently selected options.
  std::stringstream str;
  make_str_from_program_options(str);
  
  if (STRATEGY > -1)
  {
    // Save convergence graphs.   
    std::stringstream fdest, fcest, fdex, fcex, fdoutfl, fcoutfl;
    fdest << "conv_dof_est_" << str.str() << ".dat";
    fcest << "conv_cpu_est_" << str.str() << ".dat";
    fdex << "conv_dof_ex_" << str.str() << ".dat";
    fcex << "conv_cpu_ex_" << str.str() << ".dat";
    fdoutfl << "conv_dof_outfl_" << str.str() << ".dat";
    fcoutfl << "conv_cpu_outfl_" << str.str() << ".dat";
    
    graph_dof_est.save(fdest.str().c_str());
    graph_cpu_est.save(fcest.str().c_str());
    graph_dof_ex.save(fdex.str().c_str());
    graph_cpu_ex.save(fcex.str().c_str());
    graph_dof_outfl.save(fdoutfl.str().c_str());
    graph_cpu_outfl.save(fcoutfl.str().c_str());
  }
  
  if (SAVE_FINAL_SOLUTION)
  {
    // Save the final solution at 101x101 points of the domain for comparison with the semi-analytic solution.
    double y = 0.0, step = 0.01;
    int npts = int(1./step+0.5);
    double *res = new double [npts*npts];
    double *p = res;
    std::stringstream sssln;
    sssln << "sln_" << str.str() << ".dat";
    std::ofstream fs(sssln.str().c_str());
    info("Saving final solution to %s", sssln.str().c_str());
    for (int i = 0; i < npts; i++, y+=step)
    {
      std::cout << "."; std::cout.flush(); 
      double x = 0.0;
      for (int j = 0; j < npts; j++, x+=step)
        *p++ = ref_sln.get_pt_value(x, y);    
    }
    std::copy(res, res+npts*npts, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
    delete [] res;
  }
  
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
