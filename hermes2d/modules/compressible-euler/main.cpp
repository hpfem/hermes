#include "hermes2d.h"
#include "disc.h"
#include "compressible-euler.h"

int main(int argc, char* argv[])
{

  /*** RECEIVE DATA ***/

  // Check the number of command line arguments.
  if(argc != 2) error("Configuration file missing.");

  // Open configuration file.
  FILE* f = fopen(argv[1], "r");
  if(f == NULL) error("Cannot open file %s.", argv[1]);

  // Read the mesh filename.
  char * mesh_filename;
  if(!Get(f, &mesh_filename)) error("Could not read the mesh filename.");

  // Read the user setting of the usage of vector valued forms.
  bool use_vector_valued_forms;
  if(!Get(f, &use_vector_valued_forms)) error("Could not read the user setting of the usage of vector valued forms.");

  // Read the user setting of calculation of the time derivative approximation.
  bool calculate_time_derivative_approximation;
  if(!Get(f, &calculate_time_derivative_approximation)) error("Could not read the user setting of calculation of the time derivative approximation.");

  // Read the user setting of the usage of filters for visualization.
  bool filters_for_visual;
  if(!Get(f, &filters_for_visual)) error("Could not read the user setting of the usage of filters for visualization.");

  // Read the user setting of saving the solution views.
  int saving_frequency;
  if(!Get(f, &saving_frequency)) error("Could not read the user setting of saving the solution views.");

  // Read number of initial uniform mesh refinements.
  int init_ref_num;
  if(!Get(f, &init_ref_num)) error("Could not read number of initial mesh refinements.");

  // Read initial polynomial degree of elements in horizontal direction.
  int init_p_horizontal;
  if(!Get(f, &init_p_horizontal)) error("Could not read the initial polynomial degree in horizontal direction.");

  // Read initial polynomial degree of elements in vertical direction.
  int init_p_vertical;
  if(!Get(f, &init_p_vertical)) error("Could not read the initial polynomial degree in vertical direction.");

  // Read the marker for the inlet / outlet part of the boundary.
  int inlet_outlet_marker;
  if(!Get(f, &inlet_outlet_marker)) error("Could not read the marker for the inlet / outlet part of the boundary.");

  // Read the marker for the solid wall part of the boundary.
  int solid_wall_marker;
  if(!Get(f, &solid_wall_marker)) error("Could not read the marker for the solid wall part of the boundary.");

  // Read the CFL value.
  double CFL;
  if(!Get(f, &CFL)) error("Could not read the CFL value.");

  // Read the specific heat ratio.
  double kappa;
  if(!Get(f, &kappa)) error("Could not read the specific heat ratio.");

  // Read the time step.
  double tau;
  if(!Get(f, &tau)) error("Could not read the time step.");

  // Read the length of the time period.
  double T;
  if(!Get(f, &T)) error("Could not read the length of the time period.");

  // Read the exterior pressure value.
  double p_ext;
  if(!Get(f, &p_ext)) error("Could not read the exterior pressure value.");

  // Read the exterior density value.
  double rho_ext;
  if(!Get(f, &rho_ext)) error("Could not read the exterior density value.");

  // Read the exterior horizontal velocity value.
  double v1_ext;
  if(!Get(f, &v1_ext)) error("Could not read the exterior horizontal velocity value.");

  // Read the exterior vertical velocity value.
  double v2_ext;
  if(!Get(f, &v2_ext)) error("Could not read the exterior vertical velocity value.");

  /*** FEED THE DATA INTO THE COMPRESSIBLE EULER MODULE ***/

  // Initialize the Electrostatics class.
  CompressibleEuler C;

  // Set mesh filename.
  C.set_mesh_filename(std::string(mesh_filename));

  // Set the usage of vector valued forms.
  C.set_use_vector_valued_forms(use_vector_valued_forms);

  // Set the calculation of the time derivative approximation.
  
  C.set_calculate_time_derivative_approximation(calculate_time_derivative_approximation);

  // Set the usage of filters for visualization.
  C.use_filters_for_visual(filters_for_visual);

  // Set the saving of the solution views.
  C.save_every_th_solution_visualization(saving_frequency);
  
  // Set number of initial uniform mesh refinements.
  C.set_initial_mesh_refinement(init_ref_num);

  // Set initial polynomial degree of elements.
  C.set_initial_poly_degree(Ord2(init_p_horizontal, init_p_vertical));

  // Set the CFL value;
  C.set_CFL_value(CFL);

  // Set inlet / outlet boundary markers.
  C.set_inlet_outlet_marker(inlet_outlet_marker);

  // Set solid wall boundary markers.
  C.set_solid_wall_marker(solid_wall_marker);

  // Set \kappa.
  C.set_kappa(kappa);

  // Set the exterior state.
  C.set_state(p_ext, rho_ext, v1_ext, v2_ext);
  
  // Set the time interval length.
  C.set_time_interval_length(T);

  // Set the time step.
  C.set_time_step(tau);

  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e;
  bool success = C.calculate(Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  if (!success) error("Computation failed.");

  // Wait for the views to be closed.
  View::wait();

  return 1;
}
