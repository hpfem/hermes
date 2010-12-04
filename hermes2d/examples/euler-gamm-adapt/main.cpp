#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example solves the compressible Euler equations using a basic
//  piecewise-constant finite volume method.
//
//  Equations: Compressible Euler equations, perfect gas state equation.
//
//  Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
//  BC: Normal velocity component is zero on solid walls.
//      Subsonic state prescribed on inlet and outlet.
//
//  IC: Constant subsonic state identical to inlet. 
//
//  The following parameters can be changed:
 
// Experimental caching of vector valued (vector) forms.
#define HERMES_USE_VECTOR_VALUED_FORMS

const int P_INIT = 0;                     // Initial polynomial degree.                      
const int INIT_REF_NUM = 2;               // Number of initial uniform mesh refinements.                       
double CFL = 0.8;                         // CFL value.
double TAU = 1E-4;                        // Time step.

// Adaptivity.
const int UNREF_FREQ = 10;                // Every UNREF_FREQth time step the mesh is unrefined.
int REFINEMENT_COUNT = 0;                 // Number of mesh refinements between two unrefinements.
                                          // The mesh is not unrefined unless there has been a refinement since
                                          // last unrefinement.
const double THRESHOLD = 0.3;             // This is a quantitative parameter of the adapt(...) function and
                                          // it has different meanings for various adaptive strategies (see below).
const double THRESHOLD_VEL_X = 0.7;       // Threshold for adaptivity, where only second solution component is taken
                                          // into account.
const int STRATEGY = 1;                   // Adaptive strategy:
                                          // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                          //   error is processed. If more elements have similar errors, refine
                                          //   all to keep the mesh symmetric.
                                          // STRATEGY = 1 ... refine all elements whose error is larger
                                          //   than THRESHOLD times maximum element error.
                                          // STRATEGY = 2 ... refine all elements whose error is larger
                                          //   than THRESHOLD.
                                          // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;   // Predefined list of element refinement candidates. Possible values are
                                          // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                          // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                          // See User Documentation for details.
const int MESH_REGULARITY = -1;           // Maximum allowed level of hanging nodes:
                                          // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                          // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                          // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                          // Note that regular meshes are not supported, this is due to
                                          // their notoriously bad performance.
const double CONV_EXP = 1.0;              // Default value is 1.0. This parameter influences the selection of
                                          // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.75;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                          // fine mesh and coarse mesh solution in percent).
const double ERR_STOP_VEL_X = 0.5;        // Special stopping criterion for adaptivity, only second component of solution
                                          // taken into account.
const int NDOF_STOP = 100000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                          // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_SUPERLU, SOLVER_AMESOS, SOLVER_AZTECOO

// Equation parameters.
double P_EXT = 2.5;         // Exterior pressure (dimensionless).
double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
double V1_EXT = 1.25;       // Inlet x-velocity (dimensionless).
double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
double KAPPA = 1.4;         // Kappa.

double t = 0;

// Mesh boundary markers.
#define BDY_SOLID_WALL 1
#define BDY_INLET_OUTLET 2

// Numerical flux.
// For numerical fluxes, please see hermes2d/src/numerical_flux.h
NumericalFlux num_flux(KAPPA);

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

// Inlet/outlet boundary conditions.
double bc_density(double y)
{
  return RHO_EXT;
}

// Density * velocity in the x coordinate boundary condition.
double bc_density_vel_x(double y)
{
  return RHO_EXT * V1_EXT;
}

// Density * velocity in the y coordinate boundary condition.
double bc_density_vel_y(double y)
{
  return V2_EXT;
}

// Calculation of the pressure on the boundary.
double bc_pressure(double y)
{
  return P_EXT;
}

// Energy boundary condition.
double bc_energy(double y)
{
  double rho = bc_density(y);
  double rho_v_x = bc_density_vel_x(y);
  double rho_v_y = bc_density_vel_y(y);
  double pressure = bc_pressure(y);
  return pressure/(num_flux.kappa - 1.) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
}

// Calculates energy from other quantities.
// FIXME: this should be in the src/ directory, not here.
double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure)
{
  return pressure/(num_flux.kappa - 1.) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
}

// Calculates pressure from other quantities.
// FIXME: this should be in the src/ directory, not here.
double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return (num_flux.kappa - 1.) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2*rho));
}

// Calculates speed of sound.
// FIXME: this should be in the src/ directory, not here.
double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return std::sqrt(num_flux.kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy) / rho);
}

// Constant initial state (matching the supersonic inlet state).
double ic_density(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT;
}
double ic_density_vel_x(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT * V1_EXT;
}
double ic_density_vel_y(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT * V2_EXT;
}
double ic_energy(double x, double y, scalar& dx, scalar& dy)
{
  return calc_energy(RHO_EXT, RHO_EXT*V1_EXT, RHO_EXT*V2_EXT, P_EXT);
}

// Weak forms.
#include "forms.cpp"

// Filters.
#include "filters.cpp"

// Refinement criterion function.
int criterion(Element * e)
{
  if(e->vn[0]->x == 0.5 && e->vn[0]->y == 0)
    return 0;
  if(e->vn[1]->x == 0.5 && e->vn[1]->y == 0)
    return 0;
  if(e->vn[0]->x == 1.5 && e->vn[0]->y == 0)
    return 0;
  if(e->vn[1]->x == 1.5 && e->vn[1]->y == 0)
    return 0;
  
  return -1;
}



int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  basemesh.refine_by_criterion(criterion, 4);
  mesh.copy(&basemesh);

  // Initialize spaces with default shapesets.
  L2Space space_rho(&mesh, bc_types, NULL, P_INIT);
  L2Space space_rho_v_x(&mesh, bc_types, NULL, P_INIT);
  L2Space space_rho_v_y(&mesh, bc_types, NULL, P_INIT);
  L2Space space_e(&mesh, bc_types, NULL, P_INIT);

  // Initialize solutions, set initial conditions.
  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e;
  Solution rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;
  sln_rho.set_exact(&mesh, ic_density);
  sln_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  sln_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  sln_e.set_exact(&mesh, ic_energy);
  prev_rho.set_exact(&mesh, ic_density);
  prev_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  prev_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  prev_e.set_exact(&mesh, ic_energy);

  // Initialize weak formulation.
  WeakForm wf(4);

  // Bilinear forms coming from time discretization by explicit Euler's method.
  wf.add_matrix_form(0,0,callback(bilinear_form_0_0_time));
  wf.add_matrix_form(1,1,callback(bilinear_form_1_1_time));
  wf.add_matrix_form(2,2,callback(bilinear_form_2_2_time));
  wf.add_matrix_form(3,3,callback(bilinear_form_3_3_time));

  // Volumetric linear forms.
  // Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices from the previous time step.
  // First flux.
  /*
  wf.add_vector_form(0,callback(linear_form_0_1), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_x));
  wf.add_vector_form(1, callback(linear_form_1_0_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_1_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_2_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_3_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(2, callback(linear_form_2_0_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_1_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_2_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_3_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_0_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_1_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_2_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_3_first_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  // Second flux.
  
  wf.add_vector_form(0,callback(linear_form_0_2),HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_0_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_1_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_2_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(1, callback(linear_form_1_3_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(2, callback(linear_form_2_0_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_1_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_2_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
  wf.add_vector_form(2, callback(linear_form_2_3_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_0_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_1_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_2_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, callback(linear_form_3_3_second_flux), HERMES_ANY, 
                     Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  */

  // Volumetric linear forms coming from the time discretization.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form(0, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(1, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(2, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form(0,linear_form, linear_form_order, HERMES_ANY, &prev_rho);
  wf.add_vector_form(1,linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_x);
  wf.add_vector_form(2,linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_y);
  wf.add_vector_form(3,linear_form, linear_form_order, HERMES_ANY, &prev_e);
#endif

  // Surface linear forms - inner edges coming from the DG formulation.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Surface linear forms - inlet / outlet edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Surface linear forms - Solid wall edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Filters for visualization of pressure and the two components of velocity.
  SimpleFilter pressure(calc_pressure_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter u(calc_u_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter w(calc_w_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));

  // Initialize refinement selector.
  L2ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Disable weighting of refinement candidates.
  selector.set_error_weights(1, 1, 1);

  //VectorView vview("Velocity", new WinGeom(0, 0, 600, 300));
  //ScalarView sview("Pressure", new WinGeom(700, 0, 600, 300));

  ScalarView s1("w1", new WinGeom(0, 0, 620, 300));
  s1.fix_scale_width(80);
  ScalarView s2("w2", new WinGeom(625, 0, 600, 300));
  s2.fix_scale_width(50);
  ScalarView s3("w3", new WinGeom(0, 350, 620, 300));
  s3.fix_scale_width(80);
  ScalarView s4("w4", new WinGeom(625, 350, 600, 300));
  s4.fix_scale_width(50);

  // Iteration number.
  int iteration = 0;
  
  // For calculation of the time derivative of the norm of the solution approximation.
  // Not used yet in the adaptive version.
  double difference;
  double *difference_values = new double[Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e))];
  double *last_values = new double[Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e))];
  for(int i = 0; i < Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e)); i++)
      last_values[i] = 0.;
  
  // Output of the approximate time derivative.
  // Not used yet in the adaptive version.
  std::ofstream time_der_out("time_der");
  
  for(t = 0.0; t < 10; t += TAU)
  {
    info("---- Time step %d, time %3.5f.", iteration, t);

    iteration++;

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) {
      REFINEMENT_COUNT = 0;
      info("Global mesh derefinement.");
      mesh.unrefine_all_elements();
      space_rho.set_uniform_order(P_INIT);
      space_rho_v_x.set_uniform_order(P_INIT);
      space_rho_v_y.set_uniform_order(P_INIT);
      space_e.set_uniform_order(P_INIT);
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      // Global polynomial order increase = 0;
      int order_increase = 0;
      Hermes::Tuple<Space *>* ref_spaces = construct_refined_spaces(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), order_increase);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection::project_global(*ref_spaces, Hermes::Tuple<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
                     Hermes::Tuple<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), matrix_solver, 
                     Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      if(as > 1) {
        delete rsln_rho.get_mesh();
        delete rsln_rho_v_x.get_mesh();
        delete rsln_rho_v_y.get_mesh();
        delete rsln_e.get_mesh();
      }

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      bool is_linear = true;
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // The FE problem is in fact a FV problem.
      dp->set_fvm();
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
      dp->use_vector_valued_forms();
#endif
      dp->assemble(matrix, rhs);

      // Solve the linear system of the reference problem. If successful, obtain the solutions.
      if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                              Hermes::Tuple<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      else error ("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::Tuple<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
                     Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver, 
                     Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      bool solutions_for_adapt = true;
      // Error components.
      Hermes::Tuple<double> *error_components = new Hermes::Tuple<double>(4);
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
                                 Hermes::Tuple<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), solutions_for_adapt, 
                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, error_components) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // Determine the time step.
      double *solution_vector = new double[Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e))];
      OGProjection::project_global(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::Tuple<MeshFunction *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), solution_vector, matrix_solver, 
      Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double min_condition = 0;
      Element *e;
      for (int _id = 0, _max = mesh.get_max_element_id(); _id < _max; _id++) \
            if (((e) = mesh.get_element_fast(_id))->used) \
              if ((e)->active)
      {
        AsmList al;
        space_rho.get_element_assembly_list(e, &al);
        double rho = solution_vector[al.dof[0]];
        space_rho_v_x.get_element_assembly_list(e, &al);
        double v1 = solution_vector[al.dof[0]] / rho;
        space_rho_v_y.get_element_assembly_list(e, &al);
        double v2 = solution_vector[al.dof[0]] / rho;
        space_e.get_element_assembly_list(e, &al);
        double energy = solution_vector[al.dof[0]];
        
        double condition = e->get_area() / (std::sqrt(v1*v1 + v2*v2) + calc_sound_speed(rho, rho*v1, rho*v2, energy));
        
        if(condition < min_condition || min_condition == 0.)
          min_condition = condition;
      }
      if(TAU > min_condition)
        TAU = min_condition;
      if(TAU < min_condition * 0.9)
        TAU = min_condition;

      delete [] solution_vector;

      // Visualization.
      s1.show(&sln_rho);
      s2.show(&sln_rho_v_x);
      s3.show(&sln_rho_v_y);
      s4.show(&sln_e);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP && (*error_components)[1] * 100 < ERR_STOP_VEL_X) 
        done = true;
      else 
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::Tuple<RefinementSelectors::Selector *>(&selector, &selector, &selector, &selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space::get_num_dofs(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // We have to empty the cache of NeighborSearch class instances.
      NeighborSearch::empty_main_caches();

// If used, we need to clean the vector valued form caches.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
      DiscreteProblem::empty_form_caches();
#endif

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      for(int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
      delete dp;
    }
    while (done == false);

    // Debugging.
    /*    
    std::ofstream out("matrix");
    for(int i = 0; i < matrix->get_size(); i++)
      for(int j = 0; j < matrix->get_size(); j++)
        if(std::abs(matrix->get(i,j)) != 0)
          out << '(' << i << ',' << j << ')' << ':' << matrix->get(i,j) << std::endl;
    out.close();

    out.open("rhs");
      for(int j = 0; j < matrix->get_size(); j++)
        if(std::abs(rhs->get(j)) != 0)
          out << '(' << j << ')' << ':' << rhs->get(j) << std::endl;
    out.close();
     
    out.open("sol");
      for(int j = 0; j < matrix->get_size(); j++)
        out << '(' << j << ')' << ':' << solver->get_solution()[j] << std::endl;
    out.close();
    */

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);

    delete rsln_rho.get_mesh(); 
    delete rsln_rho_v_x.get_mesh(); 
    delete rsln_rho_v_y.get_mesh();
    delete rsln_e.get_mesh(); 

    // Visualization.
    /*
    pressure.reinit();
    u.reinit();
    w.reinit();
    sview.show(&pressure);
    vview.show(&u,&w);
    */
  }
  
  time_der_out.close();
  return 0;
}
