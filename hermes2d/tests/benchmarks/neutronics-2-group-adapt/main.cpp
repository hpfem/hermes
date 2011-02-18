#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"

using namespace RefinementSelectors;


/** \addtogroup t_bench_neutronics_2g Benchmarks/Neutronics 2-group
 *  \{
 *  \brief This test makes sure that the benchmark "neutronics-2-group-adapt" works correctly.
 *
 *  \section s_params Parameters
 *  - P_INIT=1 for both solution components
 *  - INIT_REF_NUM=1 for both solution components
 *  - THERSHOLD=0.3
 *  - STRATEGY=1
 *  - CAND_LIST=HP_ISO
 *  - MESH_REGULARITY=-1
 *  - ERR_STOP=0.01
 *  - CONV_EXP=1.0
 *  - NDOF_STOP=40000
 *  - ERROR_WEIGHTS=default values
 *
 *  \section s_res Expected results
 *  - DOFs: 446, 2682   (for the two solution components)
 *  - Iterations: 16    (the last iteration at which ERR_STOP is fulfilled)
 *  - Error:  3.46614%  (H1 norm of error with respect to the exact solution)
 *  - Negatives: 0      (number of negative values) 
 */
 
// INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptivity control:

const int P_INIT[2] = {1, 1};                     // Initial polynomial orders for the individual solution components.
const int INIT_REF_NUM[2] = {1, 1};               // Initial uniform mesh refinement for the individual solution components.
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const bool MULTIMESH = true;                      // true = use multi-mesh, false = use single-mesh.
                                                  // Note: in the single mesh option, the meshes are forced to be geometrically
                                                  // the same but the polynomial degrees can still vary.
const double THRESHOLD_MULTI = 0.3;               // error threshold for element refinement (multi-mesh)
const double THRESHOLD_SINGLE = 0.7;              // error threshold for element refinement (single-mesh)                                         
const CandList CAND_LIST = H2D_HP_ISO;            // Predefined list of element refinement candidates. Possible values are
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
                                                  // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 60;                     // Adaptivity process stops when the number of adaptation steps grows over
                                                  // this limit.
const int ADAPTIVITY_NORM = 2;                    // Specifies the norm used by H1Adapt to calculate the error and norm.
                                                  // ADAPTIVITY_NORM = 0 ... H1 norm.
                                                  // ADAPTIVITY_NORM = 1 ... norm defined by the diagonal parts of the bilinear form.
                                                  // ADAPTIVITY_NORM = 2 ... energy norm defined by the full (non-symmetric) bilinear form.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Variables used for reporting of results
const int ERR_PLOT = 0;         // Row in the convergence graphs for exact errors .
const int ERR_EST_PLOT = 1;     // Row in the convergence graphs for error estimates.
const int GROUP_1 = 0;          // Row in the DOF evolution graph for group 1.
const int GROUP_2 = 1;          // Row in the DOF evolution graph for group 2.


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem data:

// Two-group material properties for the 4 macro regions.
const double D[4][2]  = { {1.12, 0.6},
                          {1.2, 0.5},
                          {1.35, 0.8},
                          {1.3, 0.9}	};
const double Sr[4][2] = { {0.011, 0.13},
                          {0.09, 0.15},
                          {0.035, 0.25},
                          {0.04, 0.35}	};
const double nSf[4][2]= { {0.0025, 0.15},
                          {0.0, 0.0},
                          {0.0011, 0.1},
                          {0.004, 0.25}	};
const double chi[4][2]= { {1, 0},
                          {1, 0},
                          {1, 0},
                          {1, 0} };
const double Ss[4][2][2] = { 
                             { { 0.0, 0.0 },
                               { 0.05, 0.0 }  },
                             { { 0.0, 0.0 },
                               { 0.08, 0.0 }  },
                             { { 0.0, 0.0 },
                               { 0.025, 0.0 } },
                             { { 0.0, 0.0 },
                               { 0.014, 0.0 } } 
                           };

double a = 0., b = 1., c = (a+b)/2.;

inline int get_material(double x, double y)
{
  if (x >= a && x <= c && y >= a && y <= c) return 0;
  if (x >= c && x <= b && y >= a && y <= c) return 1;
  if (x >= c && x <= b && y >= c && y <= b) return 2;
  if (x >= a && x <= c && y >= c && y <= b) return 3;
  return -1;
}

double Q1(double x, double y)
{
  int q = get_material(x,y);

  double exfl1 = exp(-4*sqr(x))*(y/2.-sqr(y/2.));
  double exfl2 = exfl1 * (1 + sqr(sin(4*M_PI*x)) * sqr(sin(4*M_PI*y))) / 10.0;

  double L = 0.5*exp(-4*sqr(x))*(1+4*(8*sqr(x)-1)*y*(y-2))*D[q][0];
  return L + Sr[q][0]*exfl1 - chi[q][0]*nSf[q][0]*exfl1 - chi[q][1]*nSf[q][1]*exfl2;
}

double Q2(double x, double y)
{
  int q = get_material(x,y);

  double yym2 = (y-2)*y;
  double pi2 = sqr(M_PI), x2 = sqr(x), pix = M_PI*x, piy = M_PI*y;
  double cy2 = sqr(cos(4*piy)),
	 sy2 = sqr(sin(4*piy)),
	 sx2 = sqr(sin(4*pix)),
	 em4x2 = exp(-4*x2);

  double exfl1 = em4x2*(y/2.-sqr(y/2.));
  double exfl2 = exfl1 * (1 + sx2 * sy2) / 10.0;

  double L = 1./20.*em4x2*D[q][1]*(
	     1+4*(8*x2-1)*yym2+16*pi2*yym2*cy2*sx2 + 0.5*sy2*(1-4*(1+4*pi2-8*x2)*yym2 +
             (4*(1+12*pi2-8*x2)*yym2-1)*cos(8*pix) - 64*pix*yym2*sin(8*pix)) + 8*M_PI*(y-1)*sx2*sin(8*piy) );
  return L + Sr[q][1]*exfl2 - Ss[q][1][0]*exfl1;
}

double g1_D(double x, double y) {
  return 0.0;
}

double g2_D(double x, double y) {
  return 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution:

static double exact_flux1(double x, double y, double& dx, double& dy)
{
  double em4x2 = exp(-4*sqr(x));
  dx =  2.0*em4x2*x*y*(y-2);
  dy = -0.5*em4x2*(y-1);
  return em4x2*(y/2.-sqr(y/2.));
}

static double exact_flux2(double x, double y, double& dx, double& dy)
{
  double em4x2 = exp(-4*sqr(x));
  dx = 0.1*em4x2*y*(y-2)*(2*x+(2*x*sqr(sin(4*M_PI*x))-M_PI*sin(8*M_PI*x))*sqr(sin(4*M_PI*y)));
  dy = 0.05*em4x2*(1-y+sqr(sin(4*M_PI*x))*(-(y-1)*sqr(sin(4*M_PI*y))-2*M_PI*y*(y-2)*sin(8*M_PI*y)));
  return em4x2*(y/2.-sqr(y/2.)) * (1 + sqr(sin(4*M_PI*x)) * sqr(sin(4*M_PI*y))) / 10.0;
}

// APPROXIMATE SOLUTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boundary conditions:

const int BDY_FLUX = 1;
const int BDY_GAMMA = 2;
const int BDY_SYMMETRY = 3;

scalar essential_bc_values_1(double x, double y)
{
  return g1_D(x, y);
}
scalar essential_bc_values_2(double x, double y)
{
  return g2_D(x, y);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Weak forms:

#include "forms.cpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for calculating errors:

double error_total(double (*efn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*),
                   double (*nfn)(MeshFunction*, RefMap*), Hermes::vector<Solution*>& slns1, Hermes::vector<Solution*>& slns2  )
{
  double error = 0.0, norm = 0.0;

  for (unsigned int i=0; i < slns1.size(); i++) {
    error += sqr(calc_abs_error(efn, slns1[i], slns2[i]));
    if (nfn) norm += sqr(calc_norm(nfn, slns2[i]));
  }

  return (nfn ? sqrt(error/norm) : sqrt(error));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Other utility functions:

// Construct a string representing the globally set adaptivity options.
void make_str_from_adapt_opts(std::stringstream& str)
{    
  switch (CAND_LIST) {
    case H2D_H_ANISO:
    case H2D_H_ISO:
      str << "h" << P_INIT[0] << "_" << P_INIT[1];
      break;
    case H2D_P_ANISO:
    case H2D_P_ISO:
      str << "p" << INIT_REF_NUM[0] << "_" << INIT_REF_NUM[1];
      break;
    default:
      str << "hp";
      break;
  }
  switch (CAND_LIST) {
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
  
  str << (MULTIMESH ? "_multi" : "_single");
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
	Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();

  int n = 0;

  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);

    for (int i = 0; i < np; i++)
      if (uval[i] < -1e-12)
        n++;
  }

  return n;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh1, mesh2;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh1);

  if (MULTIMESH) 
  {
    // Obtain meshes for the 2nd group by cloning the mesh loaded for the 1st group.
    mesh2.copy(&mesh1);
    
    // Initial uniform refinements.
    for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
    for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();
  } 
  else // Use just one mesh for both groups.
    for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_FLUX);
  bc_types.add_bc_neumann(BDY_SYMMETRY);
  bc_types.add_bc_newton(BDY_GAMMA);

  // Enter Dirichlet boudnary values.
  BCValues bc_values_1;
  bc_values_1.add_function(BDY_FLUX, essential_bc_values_1);

  BCValues bc_values_2;
  bc_values_2.add_function(BDY_FLUX, essential_bc_values_2);

  // Create H1 space with default shapesets.
  H1Space space1(&mesh1, &bc_types, &bc_values_1, P_INIT[0]);
  H1Space space2(MULTIMESH ? &mesh2 : &mesh1, &bc_types, &bc_values_2, P_INIT[1]);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(biform_0_0), HERMES_SYM);
  wf.add_matrix_form(0, 1, callback(biform_0_1));
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(1, 1, callback(biform_1_1), HERMES_SYM);
  wf.add_vector_form(0, liform_0, liform_0_ord);
  wf.add_vector_form(1, liform_1, liform_1_ord);
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), BDY_GAMMA);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), BDY_GAMMA);
  
  // Initialize coarse and reference mesh solutions and pointers to them.
  Solution sln1, sln2; 
  Solution ref_sln1, ref_sln2;
  Hermes::vector<Solution*> slns(&sln1, &sln2);
  Hermes::vector<Solution*> ref_slns(&ref_sln1, &ref_sln2);
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  //selector.set_error_weights(2.1, 0.9, sqrt(2.0));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof("Error convergence", "Degrees of freedom", "Error [%]");
  graph_dof.add_row("exact error (H1)", "b", "-", "o");
  graph_dof.add_row("est.  error (H1)", "r", "-", "s");
  graph_dof.set_log_y();
  graph_dof.show_legend();
  graph_dof.show_grid();

  SimpleGraph graph_cpu("Error convergence", "CPU time [s]", "Error [%]");
  graph_cpu.add_row("exact error (H1)", "b", "-", "o");
  graph_cpu.add_row("est.  error (H1)", "r", "-", "s");
  graph_cpu.set_log_y();
  graph_cpu.show_legend();
  graph_cpu.show_grid();
  
  PNGGraph graph_dof_evol("Evolution of NDOF", "Adaptation step", "Num. DOF");
  graph_dof_evol.add_row("group 1", "b", "-", "o");
  graph_dof_evol.add_row("group 2", "r", "-", "s");
  graph_dof_evol.set_log_y();
  graph_dof_evol.show_legend();
  graph_dof_evol.set_legend_pos("bottom right");
  graph_dof_evol.show_grid();
    
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  double error_h1 = 0.0;
  do
  {
    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space1, &space2));

    info("!---- Adaptivity step %d ---------------------------------------------", as);

    // Construct globally refined reference meshes.
    Mesh ref_mesh1, ref_mesh2;
    ref_mesh1.copy(&mesh1);
    ref_mesh1.refine_all_elements();
    if (MULTIMESH) {
      ref_mesh2.copy(&mesh2);
      ref_mesh2.refine_all_elements();
    }
    
    // Setup spaces for the reference solution.
    int order_increase = 1;
    Space *ref_space1 = space1.dup(&ref_mesh1, order_increase);
    Space *ref_space2 = space2.dup(MULTIMESH ? &ref_mesh2 : &ref_mesh1, order_increase);

    int ref_ndof = Space::get_num_dofs(Hermes::vector<Space *>(ref_space1, ref_space2));
    info("------------------ Reference solution; NDOF=%d -------------------", ref_ndof);
                            
    if (ref_ndof >= NDOF_STOP) break;

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, Hermes::vector<Space*>(ref_space1, ref_space2), is_linear);
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    dp->assemble(matrix, rhs);
    
    // Time measurement.
    cpu_time.tick();

    // Solve the linear system associated with the reference problem.
    if(solver->solve()) 
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space*>(ref_space1, ref_space2), ref_slns);
    else error ("Matrix solver failed.\n");
    
    delete dp;
    delete ref_space1;
    delete ref_space2;
    
    // Time measurement.
    cpu_time.tick();
        
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space*>(&space1, &space2), ref_slns, slns);
    
    // View the distribution of polynomial orders on the coarse meshes.
    info("flux1_dof=%d, flux2_dof=%d", Space::get_num_dofs(&space1), Space::get_num_dofs(&space2));              
    
    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt adaptivity(Hermes::vector<Space*>(&space1, &space2));
    if (ADAPTIVITY_NORM == 2) {
      adaptivity.set_error_form(0, 0, callback(biform_0_0));
      adaptivity.set_error_form(0, 1, callback(biform_0_1));
      adaptivity.set_error_form(1, 0, callback(biform_1_0));
      adaptivity.set_error_form(1, 1, callback(biform_1_1));
    } else if (ADAPTIVITY_NORM == 1) {
      adaptivity.set_error_form(0, 0, callback(biform_0_0));
      adaptivity.set_error_form(1, 1, callback(biform_1_1));
    }
    bool solutions_for_adapt = true;
    double err_est_energ_total = adaptivity.calc_err_est(slns, ref_slns) * 100;
    
    Adapt adaptivity_proj(Hermes::vector<Space*>(&space1, &space2));

    double err_est_h1_total = adaptivity_proj.calc_err_est(slns, ref_slns) * 100;

    Hermes::vector<double> err_exact_h1;
    ExactSolution ex1(&mesh1, exact_flux1), ex2(MULTIMESH ? &mesh2 : &mesh1, exact_flux2);
    Hermes::vector<Solution*> exslns(&ex1, &ex2);
    error_h1 = adaptivity_proj.calc_err_est(slns, exslns, &err_exact_h1) * 100;
    
    // Report results.
    cpu_time.tick(); 

    // Error w.r.t. the exact solution.
    DiffFilter err_distrib_1(Hermes::vector<MeshFunction*>(&ex1, &sln1));
    DiffFilter err_distrib_2(Hermes::vector<MeshFunction*>(&ex2, &sln2));

    info("Per-component error wrt. exact solution (H1 norm): %g%%, %g%%", err_exact_h1[0] * 100, err_exact_h1[1] * 100);
    info("Total error wrt. exact solution (H1 norm): %g%%", error_h1);
    info("Total error wrt. ref. solution  (H1 norm): %g%%", err_est_h1_total);
    info("Total error wrt. ref. solution  (E norm):  %g%%", err_est_energ_total);

    if (ndof > 100) {
      // Add entry to DOF convergence graphs.
      graph_dof.add_values(ERR_PLOT, ndof, error_h1);
      graph_dof.add_values(ERR_EST_PLOT, ndof, err_est_h1_total);
      // Add entry to DOF evolution graphs.
      graph_dof_evol.add_values(GROUP_1, as, Space::get_num_dofs(&space1));
      graph_dof_evol.add_values(GROUP_2, as, Space::get_num_dofs(&space2));
      // Add entry to CPU convergence graphs.
      graph_cpu.add_values(ERR_PLOT, cpu_time.accumulated(), error_h1);
      graph_cpu.add_values(ERR_EST_PLOT, cpu_time.accumulated(), err_est_h1_total);
    }

    cpu_time.tick(HERMES_SKIP);
    
    // If err_est_energ_total too large, adapt the mesh.
    if (err_est_energ_total < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity.adapt( Hermes::vector<RefinementSelectors::Selector*>(&selector,&selector), 
                               MULTIMESH ? THRESHOLD_MULTI : THRESHOLD_SINGLE, STRATEGY, MESH_REGULARITY );
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    
    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
  }
  while (done == false);
  
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  info("Number of iterations: %d", as);
  info("NDOF: %d, %d", Space::get_num_dofs(&space1), Space::get_num_dofs(&space2));

  int n_dof_1 = Space::get_num_dofs(&space1), 
      n_dof_2 = Space::get_num_dofs(&space2);
  int n_dof_1_allowed = 300, 
      n_dof_2_allowed = 2000;
      
  int n_neg = get_num_of_neg(&sln1) + get_num_of_neg(&sln2);
  int n_neg_allowed = 0;
  
  int n_iter = as;
  int n_iter_allowed = 15;

  double error = error_h1;
  double error_allowed = 3.55;
  
  printf("n_dof_actual  = %d,%d\n", n_dof_1, n_dof_2);
  printf("n_dof_allowed = %d,%d\n", n_dof_1_allowed, n_dof_2_allowed);
  printf("n_iter_actual = %d\n", n_iter);
  printf("n_iter_allowed= %d\n", n_iter_allowed);
  printf("n_neg_actual  = %d\n", n_neg);
  printf("n_neg_allowed = %d\n", n_neg_allowed);
  printf("error_actual  = %g\n", error);
  printf("error_allowed = %g\n", error_allowed);  
  
  if (   n_dof_1 <= n_dof_1_allowed && n_dof_2 <= n_dof_2_allowed 
      && n_neg <= n_neg_allowed
      && n_iter <= n_iter_allowed
      && error <= error_allowed )   {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
