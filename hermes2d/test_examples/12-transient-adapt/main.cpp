#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Views;

//  This example is derived from example P03-timedep/03-nonlinear
//  and it shows how automatic adaptivity in space can be combined with 
//  arbitrary Runge-Kutta methods in time. The example uses fixed time 
//  step size.
//
//  For a list of available R-K methods see the file hermes_common/tables.h.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt = div[lambda(u) grad u] + f.
//
//  Nonlinearity: lambda(u) = 1 + pow(u, alpha).
//
//  Domain: square (-10, 10)^2.
//  
//  BC: Nonconstant Dirichlet.
//
//  IC: Custom initial condition matching the BC.
//
//  The following parameters can be changed:

// Hermes visualization.
const bool HERMES_VISUALIZATION = false;

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 3;                       
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;                             
// Time step. 
double time_step = 0.05;                           
// Time interval length.
const double T_FINAL = 0.51;                       

// Adaptivity
// Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_FREQ = 2;                         
// 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
const int UNREF_METHOD = 3;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Newton's method
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 20;                   

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.
ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem parameters.
// Parameter for nonlinear thermal conductivity.
const double alpha = 4.0;                         
const double heat_src = 1.0;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements(0, true);
  mesh.copy(&basemesh);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof_coarse = space.get_num_dofs();

  // Previous time level solution (initialized by initial condition).
  CustomInitialCondition sln_time_prev(&mesh);

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> f(heat_src);
  CustomWeakFormPoisson wf(&lambda, &f);

  // Next time level solution.
  Solution<double> sln_time_new(&mesh);

  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(445, 0, 410, 350));
  if(HERMES_VISUALIZATION)
  {
    view.show(&sln_time_prev);
    ordview.show(&space);
  }

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, &space, &bt);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0)
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
      case 1: mesh.copy(&basemesh);
        space.set_uniform_order(P_INIT);
        break;
      case 2: mesh.unrefine_all_elements();
        space.set_uniform_order(P_INIT);
        break;
      case 3: mesh.unrefine_all_elements();
        space.adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      }

      space.assign_dofs();
      ndof_coarse = Space<double>::get_num_dofs(&space);
    }

    // Spatial adaptivity loop. Note: sln_time_prev must not be changed
    // during spatial adaptivity.
    bool done = false; int as = 1;
    do {
      Hermes::Mixins::Loggable::Static::info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Mesh::ReferenceMeshCreator ref_mesh_creator(&mesh);
      Mesh* ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator ref_space_creator(&space, ref_mesh);
      Space<double>* ref_space = ref_space_creator.create_ref_space();
      int ndof_ref = Space<double>::get_num_dofs(ref_space);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      try
      {
        runge_kutta.set_space(ref_space);
        runge_kutta.set_verbose_output(true);
        runge_kutta.set_time(current_time);
        runge_kutta.set_time_step(time_step);
        runge_kutta.set_newton_tol(NEWTON_TOL);
        runge_kutta.rk_time_step_newton(&sln_time_prev, &sln_time_new);
      }
      catch(Exceptions::Exception& e)
      {
        std::cout << e.what();
      }

      // Project the fine mesh solution onto the coarse mesh.
      Solution<double> sln_coarse;
      Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection<double> ogProjection; ogProjection.project_global(&space, &sln_time_new, &sln_coarse);

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(&space);
      double err_est_rel_total = adaptivity->calc_err_est(&sln_coarse, &sln_time_new) * 100;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_ref: %d, err_est_rel: %g%%",
        Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space<double>::get_num_dofs(&space) >= NDOF_STOP)
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      if(HERMES_VISUALIZATION)
      {
        // Visualize the solution and mesh.
        char title[100];
        sprintf(title, "Solution<double>, time %g", current_time);
        view.set_title(title);
        view.show_mesh(false);
        view.show(&sln_time_new);
        sprintf(title, "Mesh, time %g", current_time);
        ordview.set_title(title);
        ordview.show(&space);
      }

      // Clean up.
      delete adaptivity;
      if(!done) {
        delete ref_space;
        delete sln_time_new.get_mesh();
      }
    }
    while (done == false);

    sln_time_prev.copy(&sln_time_new);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  if(HERMES_VISUALIZATION)
    View::wait();
  return 0;
}
