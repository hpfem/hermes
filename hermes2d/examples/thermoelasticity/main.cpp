#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// A massive hollow conductor is heated by induction and cooled by water running inside.
// We will model this problem using linear thermoelasticity equations, where the x-displacement,
// y-displacement, and the temperature will be approximated on individual meshes equipped
// with mutually independent adaptivity mechanisms. Use MULTI = true to use multimesh,
// MULTI = false for single-mesh (all solution components on the samemesh).
//
// PDE: Linear thermoelasticity.
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere.

const bool SOLVE_ON_COARSE_MESH = false;          // If true, coarse mesh FE problem is solved in every adaptivity step.
                                                  // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT_TEMP = 1;                        // Initial polynomial degrees in temperature mesh.
const int P_INIT_DISP = 1;                        // Initial polynomial degrees for displacement meshes.
const bool MULTI = true;                          // MULTI = true  ... use multi-mesh,
                                                  // MULTI = false ... use single-mesh.
                                                  // Note: In the single mesh option, the meshes are
                                                  // forced to be geometrically the same but the
                                                  // polynomial degrees can still vary.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 3.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double HEAT_SRC = 10000.0;                        // Heat source in the material (caused by induction heating).
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;                            // Steel: E=200 GPa.
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;                       // See http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html.

//  Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_SIDES = 2;
const int BDY_TOP = 3;
const int BDY_HOLES = 4;

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh xmesh, ymesh, tmesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &xmesh); // Master mesh.

  // Initialize multimesh hp-FEM.
  ymesh.copy(&xmesh);                  // Ydisp will share master mesh with xdisp.
  tmesh.copy(&xmesh);                  // Temp will share master mesh with xdisp.

  // Initialize boundary conditions.
  BCTypes bc_types_x_y;
  bc_types_x_y.add_bc_dirichlet(BDY_BOTTOM);
  bc_types_x_y.add_bc_neumann(Hermes::vector<int>(BDY_SIDES, BDY_TOP, BDY_HOLES));

  BCTypes bc_types_t;
  bc_types_t.add_bc_dirichlet(BDY_HOLES);
  bc_types_t.add_bc_neumann(Hermes::vector<int>(BDY_SIDES, BDY_TOP, BDY_BOTTOM)); 

  // Enter Dirichlet boundary values.
  BCValues bc_values_x_y;
  bc_values_x_y.add_zero(BDY_BOTTOM);

  BCValues bc_values_t;
  bc_values_t.add_const(BDY_HOLES, TEMP_INNER);

  // Create H1 spaces with default shapesets.
  H1Space xdisp(&xmesh, &bc_types_x_y, &bc_values_x_y, P_INIT_DISP);
  H1Space ydisp(MULTI ? &ymesh : &xmesh, &bc_types_x_y, &bc_values_x_y, P_INIT_DISP);
  H1Space temp(MULTI ? &tmesh : &xmesh, &bc_types_t, &bc_values_t, P_INIT_TEMP);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), HERMES_SYM);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2));
  wf.add_vector_form(1, callback(linear_form_1));
  wf.add_vector_form(2, callback(linear_form_2));
  wf.add_vector_form_surf(2, callback(linear_form_surf_2));

  // Initialize coarse and reference mesh solutions.
  Solution xdisp_sln, ydisp_sln, temp_sln, ref_xdisp_sln, ref_ydisp_sln, ref_temp_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view_0("Solution[xdisp]", new WinGeom(0, 0, 450, 350));
  s_view_0.show_mesh(false);
  ScalarView s_view_1("Solution[ydisp]", new WinGeom(460, 0, 450, 350));
  s_view_1.show_mesh(false);
  ScalarView s_view_2("Solution[temp]", new WinGeom(920, 0, 450, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_0("Mesh[xdisp]", new WinGeom(0, 360, 450, 350));
  OrderView  o_view_1("Mesh[ydisp]", new WinGeom(460, 360, 450, 350));
  OrderView  o_view_2("Mesh[temp]", new WinGeom(920, 360, 450, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&xdisp, &ydisp, &temp));

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. If successful, obtain the solutions.
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                            Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln));
    else error ("Matrix solver failed.\n");
  
    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space *>(&xdisp, &ydisp, &temp), Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln), 
                   Hermes::vector<Solution *>(&xdisp_sln, &ydisp_sln, &temp_sln), matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&xdisp_sln); 
    o_view_0.show(&xdisp);
    s_view_1.show(&ydisp_sln); 
    o_view_1.show(&ydisp);
    s_view_2.show(&temp_sln); 
    o_view_2.show(&temp);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&xdisp, &ydisp, &temp));
    adaptivity->set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    adaptivity->set_error_form(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    adaptivity->set_error_form(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
    adaptivity->set_error_form(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    adaptivity->set_error_form(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    adaptivity->set_error_form(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
    adaptivity->set_error_form(2, 2, bilinear_form_2_2<scalar, scalar>, bilinear_form_2_2<Ord, Ord>);

    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&xdisp_sln, &ydisp_sln, &temp_sln), 
                              Hermes::vector<Solution *>(&ref_xdisp_sln, &ref_ydisp_sln, &ref_temp_sln), &err_est_rel) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[xdisp]: %d, ndof_fine[xdisp]: %d, err_est_rel[xdisp]: %g%%", 
         xdisp.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[0]), err_est_rel[0]*100);
    info("ndof_coarse[ydisp]: %d, ndof_fine[ydisp]: %d, err_est_rel[ydisp]: %g%%",
         ydisp.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[1]), err_est_rel[1]*100);
    info("ndof_coarse[temp]: %d, ndof_fine[temp]: %d, err_est_rel[temp]: %g%%",
         temp.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[2]), err_est_rel[2]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel_total: %g%%",
         Space::get_num_dofs(Hermes::vector<Space *>(&xdisp, &ydisp, &temp)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(Hermes::vector<Space *>(&xdisp, &ydisp, &temp)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector), 
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space::get_num_dofs(Hermes::vector<Space *>(&xdisp, &ydisp, &temp)) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
    delete dp;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  s_view_0.set_title("Fine mesh Solution[xdisp]");
  s_view_0.show(&ref_xdisp_sln);
  s_view_1.set_title("Fine mesh Solution[ydisp]");
  s_view_1.show(&ref_ydisp_sln);
  s_view_1.set_title("Fine mesh Solution[temp]");
  s_view_1.show(&ref_temp_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
};

