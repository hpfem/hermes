#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
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
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

const int P_INIT_XVEL = 2;
const int P_INIT_YVEL = 2;
const int P_INIT_PRESS = 1;
const int P_INIT_LSET = 1;

// Initial conditions for velocity, pressure, and level-set function.
BCType xvel_bc_types(int marker)
{ 
/*  if (marker == 1) return BC_ESSENTIAL; 
  else return BC_NONE;*/
    return BC_ESSENTIAL;
}


BCType yvel_bc_types(int marker)
{ 
/*  if ((marker == 2) || (marker == 3)) return BC_ESSENTIAL; 
  else return BC_NONE;*/
     return BC_ESSENTIAL;
}


BCType lset_bc_types(int marker)
{ 
//  return BC_NONE; 
  return BC_ESSENTIAL;
}

scalar essential_bc_values_lset(int essential_marker, double x, double y)
{
/*  if ((marker == 1) && (x < 0.5)) return 0.5;
  if ((marker == 1) && (x > 0.5)) return -0.5;
  if ((marker == 2) || (marker == 3)) return 0.5 - x;*/
  return 0.0;
}


BCType press_bc_types(int marker)
{ 
  return BC_NONE;
}

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
  double elem_error[ne];
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
    mesh->refine_element(elist[i].id);
    for (j = 0; j < 4; j++)
      space->set_element_order(e->sons[j]->id, space->get_element_order(elist[i].id));

    processed_error += elist[i].error;  
  }    

  delete [] elist;
 
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/*
bool adapt_solution(H1OrthoHP* h1hp, Mesh* mesh, Mesh* rmesh, Solution* sln, Solution* rsln, Space* space)
{
  int i, j;
  
  int ne = mesh->get_max_element_id() + 1;
  double elem_error[ne];
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
 
  // sort elements by their error
  qsort(elist, n, sizeof(ElemList), compare_error);
  
   // refine worst elements
  double processed_error = 0.0;
  double err = 1000.0;

  for (i = 0; i < n; i++)
  {
    if ((processed_error > 0.5 * (total_err)) && (fabs((elist[i].error - err)/err) > 1e-3))
      break;
    err = elist[i].error;

    e = mesh->get_element(elist[i].id);

    int split = 0;
    int p[4];
    int current = space->get_element_order(elist[i].id);
    rsln->set_quad_2d(&g_quad_2d_std);   
    rsln->enable_transform(false);
    h1hp->get_optimal_refinement(e, current, rsln, split, p);
    rsln->enable_transform(true);

    if (split == -1)
      space->set_element_order(elist[i].id, p[0]);
    else if (split == 0){
      mesh->refine_element(elist[i].id);
      for (j = 0; j < 4; j++)
        space->set_element_order(e->sons[j]->id, p[j]);
    }
    else {
      mesh->refine_element(elist[i].id, split);
      for (j = 0; j < 2; j++)
        space->set_element_order(e->sons[ (split == 1) ? j : j+2 ]->id, p[j]);
    }    

    processed_error += elist[i].error;  
  }    

  delete [] elist;
 
  return false;
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////

/*void refine_interface(Mesh* mesh, Solution* l)
{
  Element* e;
  int ne = mesh->get_max_element_id() + 1;
  int list[ne];
  int k = 0;
  for_all_active_elements(e, mesh)
  {
    l->set_active_element(e);
    Quad2D* quad = l->get_quad_2d();
    int o = 10;
    l->set_quad_order(o, FN_VAL);

    scalar* lval = l->get_fn_values();
    int np = quad->get_num_points(o);
    bool found = false;
    bool sign = (lval[0] < 0.0) ? false : true;
    for (int i = 1; i < np; i++)
    {
      bool lsign = (lval[i] < 0.0) ? false : true;
      if (lsign != sign) { found = true; break; }
    }
    if (found)
      list[k++] = e->id;
  }
  for (int i = 0; i < k; i++)
    mesh->refine_element(list[i]);
  

}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////
/*
static void init_mesh(Mesh* mesh, Space* space)
{
  Mesh old_mesh;
  old_mesh.copy(mesh);
  mesh->unrefine_all_elements(); 

  int ne = old_mesh.get_max_element_id() + 1;
  int orders[ne];  

  Element* e;
  for_all_active_elements(e, mesh)
  {
    Element* p = old_mesh.get_element(e->id);
    int o;
    if (p->active)
    {
      o = get_h_order(space->get_element_order(p->id));
    }
    else
    {
      o = 10;
      for (int i = 0; i < 4; i++)
        if (p->sons[i] != NULL)
          o = std::min(o, get_h_order(space->get_element_order(p->sons[i]->id)));
    }
    orders[e->id] = o;
  }
  for_all_active_elements(e, mesh)
  {
    space->set_element_order(e->id, std::max(1, orders[e->id]-1));
  }
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////


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

  // Spaces for velocities and pressure.
  H1Space xvel(&mesh1, xvel_bc_types, NULL, P_INIT_XVEL);
  H1Space yvel(&mesh1, yvel_bc_types, NULL, P_INIT_YVEL);
  H1Space press(&mesh1, press_bc_types, NULL, P_INIT_PRESS);
  H1Space lset(&mesh2, xvel_bc_types, essential_bc_values_lset, P_INIT_LSET);
  int ndof_1 = xvel.get_num_dofs();
  int ndof_2 = yvel.get_num_dofs();
  int ndof_3 = press.get_num_dofs();
  int ndof_4 = lset.get_num_dofs();
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
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&lprev, &xprev, &yprev));

  // OLD: wf.add_biform(0, 2, bilinear_form_0_2);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2), H2D_UNSYM, H2D_ANY);

  // OLD: wf.add_biform(1, 1, bilinear_form_0_0, 0, 0, 3, &lprev, &xprev, &yprev);
  wf.add_matrix_form(1, 1, callback(bilinear_form_0_0), H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&lprev, &xprev, &yprev));

  // OLD: wf.add_biform(1, 2, bilinear_form_1_2);
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2), H2D_UNSYM, H2D_ANY);

  // OLD: wf.add_biform(2, 0, bilinear_form_2_0);
  wf.add_matrix_form(2, 0, callback(bilinear_form_2_0), H2D_UNSYM, H2D_ANY);

  // OLD: wf.add_biform(2, 1, bilinear_form_2_1);
  wf.add_matrix_form(2, 1, callback(bilinear_form_2_1), H2D_UNSYM, H2D_ANY);

  // OLD: wf.add_biform(3, 3, bilinear_form_3_3, 0, 0, 2, &xprev, &yprev);
  wf.add_matrix_form(3, 3, callback(bilinear_form_3_3), H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xprev, &yprev));

  // OLD: wf.add_liform(0, linear_form_0, 0, 2, &lprev, &xprev);
  wf.add_vector_form(0, callback(linear_form_0), H2D_ANY, Tuple<MeshFunction*>(&lprev, &yprev));

  // OLD: wf.add_liform(1, linear_form_1, 0, 2, &lprev, &yprev);
  wf.add_vector_form(1, callback(linear_form_1), H2D_ANY, Tuple<MeshFunction*>(&lprev, &yprev));

  // OLD: wf.add_liform(3, linear_form_3, 0, 1, &lprev);
  wf.add_vector_form(3, callback(linear_form_3), H2D_ANY, &lprev);

  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  Solution u1, u2, u3, u4;
  Solution r1, r2, r3, r4;

  // Initialize adaptivity parameters.
  double to_be_processed = 0;
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY,
                          MESH_REGULARITY, to_be_processed, H2D_TOTAL_ERROR_REL, H2D_ELEMENT_ERROR_REL);
  apt.set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
  apt.set_error_form(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
  apt.set_error_form(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
  apt.set_error_form(2, 0, bilinear_form_2_0<scalar, scalar>, bilinear_form_2_0<Ord, Ord>);
  apt.set_error_form(2, 1, bilinear_form_2_1<scalar, scalar>, bilinear_form_2_1<Ord, Ord>);
  apt.set_error_form(3, 3, bilinear_form_3_3<scalar, scalar>, bilinear_form_3_3<Ord, Ord>);

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

  bool verbose = true; 
  solve_linear_adapt(Tuple<Space *>(&xvel, &yvel, &press, &lset), &wf, NULL, matrix_solver,
                     Tuple<int>(H2D_H1_NORM, H2D_H1_NORM, H2D_H1_NORM, H2D_H1_NORM),
                     Tuple<Solution *>(&u1, &u2, &u3, &u4),
                     Tuple<Solution *>(&r1, &r2, &r3, &r4),
                     Tuple<WinGeom *>(), Tuple<WinGeom *>(),// Do not show solutions or meshes.
                     Tuple<RefinementSelectors::Selector *> (&selector, &selector, &selector, &selector), 
                     &apt, verbose);

  u1_sln_view.show(&u1);
  u2_sln_view.show(&u2);
  u3_sln_view.show(&u3);
  u4_sln_view.show(&u4);
  u1_order_view.show(&xvel);
  u2_order_view.show(&yvel);
  u3_order_view.show(&press);
  u4_order_view.show(&lset);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
