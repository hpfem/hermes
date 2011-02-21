#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

const double tau = 0.05;   // Time step.
//const double Ro2 = 1000.0; // The density, water.
const double Ro1 = 1000.0; // The density.
//const double Ro1 = 1060.0;  // The density, blood.
//const double Ro2 = 13.534; // The density, mercury.
const double Ro2 = 900.0;  // The density, oil.
//const double Nu2 = 0.001; // The viscosity, water.
const double Nu1 = 0.01; // The viscosity.
//const double Nu2 = 0.00153; // The viscosity, mercury.
const double Nu2 = 0.081; // The viscosity, oil.
//const double Nu1 = 0.003; // The viscosity, blood.

const CandList CAND_LIST = H2D_HP_ANISO;
const int MESH_REGULARITY = -1;
const double ERR_STOP = 0.5;
const int NDOF_STOP = 60000; 
const double CONV_EXP = 1.0;
const double THRESHOLD = 0.3;
const bool MULTI = true; 
const int STRATEGY = 1; 
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const int P_INIT_XVEL = 2;
const int P_INIT_YVEL = 2;
const int P_INIT_PRESS = 1;
const int P_INIT_LSET = 1;

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Weak forms.
#include "forms.cpp"

scalar exactfn(double x, double y, scalar& dx , scalar& dy) { 
/*  dx=-1; dy=0; //return (y<.48) ? -1 : 1;
  return -x + 0.5;*/
/*  dx = 4.0 - 8*x; dy = 4.0 - 8*y;
  return 1.0 - 5*(sqr(2*x - 1) + sqr(2*y - 1)) ;*/
  dx = 0; dy = 0;
  return sin(2.0 * M_PI * x) * sin(M_PI * y);

}

// Velocities from the previous time step.
Solution xprev, yprev, lprev;

struct ElemList
{
  int id;
  double error; /// Element error.
};

// Helper function for sorting elements by their error.
static int compare_error(const void* x1, const void* x2)
{
  const ElemList *ei1 = (const ElemList*) x1;
  const ElemList *ei2 = (const ElemList*) x2;
  return (ei1)->error < (ei2)->error ? 1 : -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool adapt_velocity(Mesh* mesh,  Mesh* rmesh, MeshFunction* sln, MeshFunction* rsln, Space* space)
{
  int i, j;
  
  int ne = mesh->get_max_element_id() + 1;
  double* elem_error = new double[ne];
  memset(elem_error, 0, sizeof(double) * ne);
  
  double error;

  int num = mesh->get_num_active_elements();
  ElemList* elist = new ElemList[num]; 

  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  rsln->set_quad_2d(quad);
    
  Mesh* meshes[2] = { mesh, rmesh };
  Transformable* tr[2] = { sln, rsln };
  Traverse trav;
  trav.begin(2, meshes, tr);

  Element** ee;
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    RefMap* ru = sln->get_refmap();
    RefMap* rr = rsln->get_refmap();
    elem_error[ee[0]->id] += int_l2_error(sln, rsln, ru, rr);
    
  }
  trav.finish();

  Element* e;     
  double total_err = 0.0;
  int n = 0;
  for_all_active_elements(e, mesh)
  {
    elist[n].id = e->id;
    elist[n].error = elem_error[e->id];
    total_err += elem_error[e->id];
    n++;
  }
 
  // Sort elements by their error.
  qsort(elist, n, sizeof(ElemList), compare_error);
  
   // Refine worst elements.
  double processed_error = 0.0;
  double err = 1000.0;

  for (i = 0; i < n; i++)
  {
    if ((processed_error > 0.5 * (total_err)) && (fabs((elist[i].error - err)/err) > 1e-3))
      break;
    err = elist[i].error;

    e = mesh->get_element(elist[i].id);
    mesh->refine_element_id(elist[i].id);
    for (j = 0; j < 4; j++)
      space->set_element_order(e->sons[j]->id, space->get_element_order(elist[i].id));

    processed_error += elist[i].error;  
  }    

  delete [] elist;
  delete [] elem_error;
 
  return false;
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh file.
  Mesh mesh1;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh1);
  mesh1.refine_all_elements();
  mesh1.refine_all_elements();
  mesh1.refine_all_elements();

  Mesh mesh2;
  mloader.load("square_quad.mesh", &mesh2);
  mesh2.refine_all_elements();
  mesh2.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  BCTypes bc_types_press;
  bc_types_press.add_bc_none(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Spaces for velocities and pressure.
  H1Space xvel(&mesh1, &bc_types, &bc_values, P_INIT_XVEL);
  H1Space yvel(&mesh1, &bc_types, &bc_values, P_INIT_YVEL);
  H1Space press(&mesh1, &bc_types_press, P_INIT_PRESS);
  H1Space lset(&mesh2, &bc_types, &bc_values, P_INIT_LSET);

  int ndof_1 = xvel.Space::get_num_dofs();
  int ndof_2 = yvel.Space::get_num_dofs();
  int ndof_3 = press.Space::get_num_dofs();
  int ndof_4 = lset.Space::get_num_dofs();
  info("ndof = %d (xvel), %d (yvel), %d (press), %d (lset)", ndof_1, ndof_2, ndof_3, ndof_4);

  // Zero initial xprev and yprev values.
  xprev.set_zero(&mesh1);
  yprev.set_zero(&mesh1);
  lprev.set_exact(&mesh2, exactfn); 

  // OLD: Reference solution
  //  Mesh rmesh1, rmesh2;

  // OLD: 
  //  H1Space rspace1(&rmesh1, xvel_bc_types, NULL, P_INIT_XVEL);
  //  H1Space rspace2(&rmesh1, yvel_bc_types, NULL, P_INIT_YVEL);
  //  H1Space rspace3(&rmesh1, press_bc_types, NULL, P_INIT_PRESS);
  //  H1Space rspace4(&rmesh2, lset_bc_types, essential_bc_values_lset, P_INIT_LSET);

  // Problem initalization
  WeakForm wf(4);
  // OLD: wf.add_biform(0, 0, bilinear_form_0_0, 0, 0, 3, &lprev, &xprev, &yprev);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), HERMES_NONSYM, HERMES_ANY, Hermes::vector<MeshFunction*>(&lprev, &xprev, &yprev));

  // OLD: wf.add_biform(0, 2, bilinear_form_0_2);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2), HERMES_NONSYM, HERMES_ANY);

  // OLD: wf.add_biform(1, 1, bilinear_form_0_0, 0, 0, 3, &lprev, &xprev, &yprev);
  wf.add_matrix_form(1, 1, callback(bilinear_form_0_0), HERMES_NONSYM, HERMES_ANY, Hermes::vector<MeshFunction*>(&lprev, &xprev, &yprev));

  // OLD: wf.add_biform(1, 2, bilinear_form_1_2);
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2), HERMES_NONSYM, HERMES_ANY);

  // OLD: wf.add_biform(2, 0, bilinear_form_2_0);
  wf.add_matrix_form(2, 0, callback(bilinear_form_2_0), HERMES_NONSYM, HERMES_ANY);

  // OLD: wf.add_biform(2, 1, bilinear_form_2_1);
  wf.add_matrix_form(2, 1, callback(bilinear_form_2_1), HERMES_NONSYM, HERMES_ANY);

  // OLD: wf.add_biform(3, 3, bilinear_form_3_3, 0, 0, 2, &xprev, &yprev);
  wf.add_matrix_form(3, 3, callback(bilinear_form_3_3), HERMES_NONSYM, HERMES_ANY, Hermes::vector<MeshFunction*>(&xprev, &yprev));

  // OLD: wf.add_liform(0, linear_form_0, 0, 2, &lprev, &xprev);
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, Hermes::vector<MeshFunction*>(&lprev, &yprev));

  // OLD: wf.add_liform(1, linear_form_1, 0, 2, &lprev, &yprev);
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, Hermes::vector<MeshFunction*>(&lprev, &yprev));

  // OLD: wf.add_liform(3, linear_form_3, 0, 1, &lprev);
  wf.add_vector_form(3, callback(linear_form_3), HERMES_ANY, &lprev);

  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  Solution u1, u2, u3, u4;
  Solution r1, r2, r3, r4;

  // Initialize adaptivity parameters.
  double to_be_processed = 0;
  

  // Geometry and position of visualization windows.
  WinGeom* u1_sln_win_geom = new WinGeom(0, 0, 300, 450);
  WinGeom* u2_sln_win_geom = new WinGeom(0, 0, 300, 450);
  WinGeom* u3_sln_win_geom = new WinGeom(310, 0, 300, 450);
  WinGeom* u4_sln_win_geom = new WinGeom(310, 0, 300, 450);
  WinGeom* u1_mesh_win_geom = new WinGeom(620, 0, 280, 450);
  WinGeom* u2_mesh_win_geom = new WinGeom(620, 0, 280, 450);
  WinGeom* u3_mesh_win_geom = new WinGeom(910, 0, 280, 450);
  WinGeom* u4_mesh_win_geom = new WinGeom(910, 0, 280, 450);

  // Initialize views.
  ScalarView u1_sln_view("[1]", u1_sln_win_geom);
  ScalarView u2_sln_view("[2]", u2_sln_win_geom);
  ScalarView u3_sln_view("[3]", u3_sln_win_geom);
  ScalarView u4_sln_view("[4]", u4_sln_win_geom);
  OrderView u1_order_view("[1]", u1_mesh_win_geom);
  OrderView u2_order_view("[2]", u2_mesh_win_geom);
  OrderView u3_order_view("[3]", u3_mesh_win_geom);
  OrderView u4_order_view("[4]", u4_mesh_win_geom);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  cpu_time.tick();

  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset));

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem dp(&wf, *ref_spaces, is_linear);
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    dp.assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, Hermes::vector<Solution *>(&r1, &r2, &r3, &r4));
    else error ("Matrix solver failed.\n");

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset), Hermes::vector<Solution *>(&r1, &r2, &r3, &r4), Hermes::vector<Solution *>(&u1, &u2, &u3, &u4), matrix_solver); 

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution(s).
    u1_sln_view.show(&u1);
    u2_sln_view.show(&u2);
    u3_sln_view.show(&u3);
    u4_sln_view.show(&u4);
    u1_order_view.show(&xvel);
    u2_order_view.show(&yvel);
    u3_order_view.show(&press);
    u4_order_view.show(&lset);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset), Hermes::vector<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM, HERMES_H1_NORM, HERMES_H1_NORM));
    /*
    adaptivity->set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    adaptivity->set_error_form(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
    adaptivity->set_error_form(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
    adaptivity->set_error_form(2, 0, bilinear_form_2_0<scalar, scalar>, bilinear_form_2_0<Ord, Ord>);
    adaptivity->set_error_form(2, 1, bilinear_form_2_1<scalar, scalar>, bilinear_form_2_1<Ord, Ord>);
    adaptivity->set_error_form(3, 3, bilinear_form_3_3<scalar, scalar>, bilinear_form_3_3<Ord, Ord>);
    */
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&u1, &u2, &u3, &u4), 
                               Hermes::vector<Solution *>(&r1, &r2, &r3, &r4), &err_est_rel) * 100;

    // Report results.
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
         Space::get_num_dofs(&xvel), Space::get_num_dofs((*ref_spaces)[0]));
    info("err_est_rel[0]: %g%%", err_est_rel[0]*100);

    info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
         Space::get_num_dofs(&yvel), Space::get_num_dofs((*ref_spaces)[1]));
    info("err_est_rel[1]: %g%%", err_est_rel[1]*100);

    info("ndof_coarse[2]: %d, ndof_fine[2]: %d",
         Space::get_num_dofs(&press), Space::get_num_dofs((*ref_spaces)[2]));
    info("err_est_rel[2]: %g%%", err_est_rel[2]*100);

    info("ndof_coarse[3]: %d, ndof_fine[3]: %d",
         Space::get_num_dofs(&lset), Space::get_num_dofs((*ref_spaces)[3]));
    info("err_est_rel[3]: %g%%", err_est_rel[3]*100);


    info("ndof_coarse_total: %d, ndof_fine_total: %d",
         Space::get_num_dofs(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset)), Space::get_num_dofs(*ref_spaces));
    info("err_est_rel_total: %g%%", err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *> (&selector, &selector, &selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY, to_be_processed);

      if (Space::get_num_dofs(Hermes::vector<Space *>(&xvel, &yvel, &press, &lset)) >= NDOF_STOP) done = true;
    }

    delete matrix;
    delete rhs;
    delete solver;
    delete adaptivity;
    as++;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
