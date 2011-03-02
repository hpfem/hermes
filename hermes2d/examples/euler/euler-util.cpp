extern NumericalFlux num_flux;

// Calculates energy from other quantities.
double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure)
{
  return pressure/(num_flux.kappa - 1.) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
};

// Calculates pressure from other quantities.
double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return (num_flux.kappa - 1.) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2*rho));
};

// Calculates speed of sound.
double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return std::sqrt(num_flux.kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy) / rho);
};

class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(Hermes::vector<Space *> spaces, 
                        Hermes::vector<Solution *> solutions);

  /// Destructor.
   ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids(double threshold);

  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  double calculate_jumps(Element* e, int edge_i);

  /// Calculates the norm of the solution on the central element.
  double calculate_norm(Element* e, int edge_i);

protected:
  /// Members.
  Hermes::vector<Space *> spaces;
  Hermes::vector<Solution *> solutions;
  Mesh* mesh;
  std::set<int> discontinuous_element_ids;
};


DiscontinuityDetector::DiscontinuityDetector(Hermes::vector<Space *> spaces, 
                        Hermes::vector<Solution *> solutions) : spaces(spaces), solutions(solutions)
{
  // A check that all meshes are the same in the spaces.
  Mesh* mesh0 = spaces[0]->get_mesh();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh() != mesh0)
      error("So far DiscontinuityDetector works only for single mesh.");
  mesh = mesh0;
};

DiscontinuityDetector::~DiscontinuityDetector()
{};

std::set<int>& DiscontinuityDetector::get_discontinuous_element_ids(double threshold)
{
  Element* e;
  for_all_active_elements(e, mesh) {
    for(int edge_i = 0; edge_i < e->get_num_surf(); edge_i++)
      if(calculate_relative_flow_direction(e, edge_i) < -1e-3 && !e->en[edge_i]->bnd) {
        double jump = calculate_jumps(e, edge_i);
        double diameter_indicator = std::pow(e->get_diameter(), 
                                             (0.5 * (H2D_GET_H_ORDER(spaces[0]->get_element_order(e->id)) 
                                                     + 
                                                     H2D_GET_V_ORDER(spaces[0]->get_element_order(e->id)))
                                              + 1)
                                              / 2);
        double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % 4]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % 4]->y - e->vn[edge_i]->y, 2));
        double norm = calculate_norm(e, edge_i);
        double discontinuity_detector = jump / (diameter_indicator * edge_length * norm);
        if(discontinuity_detector > threshold)
          discontinuous_element_ids.insert(e->id);
      }
  }
  return discontinuous_element_ids;
};

double DiscontinuityDetector::calculate_relative_flow_direction(Element* e, int edge_i)
{
  // Set active element to the two solutions (density_vel_x, density_vel_y).
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 5);
  double3* pt = solutions[1]->get_quad_2d()->get_points(eo);
  int np = solutions[1]->get_quad_2d()->get_num_points(eo);

  Geom<double>* geom = init_geom_surf(solutions[1]->get_refmap(), &surf_pos, eo);
  double3* tan = solutions[1]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
      jwt[i] = pt[i][2] * tan[i][2];
 
  // Calculate.
  Func<scalar>* density_vel_x = init_fn(solutions[1], eo);
  Func<scalar>* density_vel_y = init_fn(solutions[2], eo);

  double result = 0.0;
  for(int point_i = 0; point_i < np; point_i++)
    result += jwt[point_i] * density_vel_x->val[point_i] * geom->nx[point_i] + density_vel_y->val[point_i] * geom->ny[point_i];

  return result;
};

double DiscontinuityDetector::calculate_jumps(Element* e, int edge_i)
{
  

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 5);
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo);
  int np = solutions[0]->get_quad_2d()->get_num_points(eo);

  // Initialize the NeighborSearch.
  NeighborSearch ns(e, mesh);
  ns.set_active_edge(edge_i);

  // The value to be returned.
  double result = 0.0;

  // Go through all neighbors.
  for(unsigned int neighbor_i = 0; neighbor_i < ns.get_num_neighbors(); neighbor_i++) {
    ns.active_segment = neighbor_i;
    ns.neighb_el = ns.neighbors[neighbor_i];
    ns.neighbor_edge = ns.neighbor_edges[neighbor_i];

    // Set active element to the solutions.
    solutions[0]->set_active_element(e);
    solutions[1]->set_active_element(e);
    solutions[2]->set_active_element(e);
    solutions[3]->set_active_element(e);
  
    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.central_n_trans[neighbor_i]; trf_i++) {
      solutions[0]->push_transform(ns.central_transformations[neighbor_i][trf_i]);
      solutions[1]->push_transform(ns.central_transformations[neighbor_i][trf_i]);
      solutions[2]->push_transform(ns.central_transformations[neighbor_i][trf_i]);
      solutions[3]->push_transform(ns.central_transformations[neighbor_i][trf_i]);
    }

    Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), &surf_pos, eo);
    double3* tan = solutions[0]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
    double* jwt = new double[np];
    for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * tan[i][2];
 
    // Prepare functions on the central element.
    Func<scalar>* density = init_fn(solutions[0], eo);
    Func<scalar>* density_vel_x = init_fn(solutions[1], eo);
    Func<scalar>* density_vel_y = init_fn(solutions[2], eo);
    Func<scalar>* density_energy = init_fn(solutions[3], eo);

    // Set neighbor element to the solutions.
    solutions[0]->set_active_element(ns.neighb_el);
    solutions[1]->set_active_element(ns.neighb_el);
    solutions[2]->set_active_element(ns.neighb_el);
    solutions[3]->set_active_element(ns.neighb_el);

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.neighbor_n_trans[neighbor_i]; trf_i++) {
      solutions[0]->push_transform(ns.neighbor_transformations[neighbor_i][trf_i]);
      solutions[1]->push_transform(ns.neighbor_transformations[neighbor_i][trf_i]);
      solutions[2]->push_transform(ns.neighbor_transformations[neighbor_i][trf_i]);
      solutions[3]->push_transform(ns.neighbor_transformations[neighbor_i][trf_i]);
    }

    // Prepare functions on the neighbor element.
    Func<scalar>* density_neighbor = init_fn(solutions[0], eo);
    Func<scalar>* density_vel_x_neighbor = init_fn(solutions[1], eo);
    Func<scalar>* density_vel_y_neighbor = init_fn(solutions[2], eo);
    Func<scalar>* density_energy_neighbor = init_fn(solutions[3], eo);

    DiscontinuousFunc<scalar> density_discontinuous(density, density_neighbor, true);
    DiscontinuousFunc<scalar> density_vel_x_discontinuous(density_vel_x, density_vel_x_neighbor, true);
    DiscontinuousFunc<scalar> density_vel_y_discontinuous(density_vel_y, density_vel_y_neighbor, true);
    DiscontinuousFunc<scalar> density_energy_discontinuous(density_energy, density_energy_neighbor, true);

    for(int point_i = 0; point_i < np; point_i++)
      result += jwt[point_i] * (
      std::pow(density_discontinuous.get_val_central(point_i) - density_discontinuous.get_val_neighbor(point_i), 2) + 
      std::pow(density_vel_x_discontinuous.get_val_central(point_i) - density_vel_x_discontinuous.get_val_neighbor(point_i), 2) + 
      std::pow(density_vel_y_discontinuous.get_val_central(point_i) - density_vel_y_discontinuous.get_val_neighbor(point_i), 2) + 
      std::pow(density_energy_discontinuous.get_val_central(point_i) - density_energy_discontinuous.get_val_neighbor(point_i), 2));
  }

  return std::sqrt(result);
};

double DiscontinuityDetector::calculate_norm(Element* e, int edge_i)
{
  // Set active element to the solutions.
  solutions[0]->set_active_element(e);
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);
  solutions[3]->set_active_element(e);

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 5);
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo);
  int np = solutions[0]->get_quad_2d()->get_num_points(eo);

  Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), &surf_pos, eo);
  double3* tan = solutions[0]->get_refmap()->get_tangent(surf_pos.surf_num, eo);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
      jwt[i] = pt[i][2] * tan[i][2];
 
  // Calculate.
  Func<scalar>* density = init_fn(solutions[0], eo);
  Func<scalar>* density_vel_x = init_fn(solutions[1], eo);
  Func<scalar>* density_vel_y = init_fn(solutions[2], eo);
  Func<scalar>* density_energy = init_fn(solutions[3], eo);

  double result = 0.0;
  for(int point_i = 0; point_i < np; point_i++)
    result += jwt[point_i] * density->val[point_i] * density->val[point_i] + density_vel_x->val[point_i] * density_vel_x->val[point_i]
              + density_vel_y->val[point_i] * density_vel_y->val[point_i] + density_energy->val[point_i] * density_energy->val[point_i];

  return std::sqrt(result);

};

class FluxLimiter
{
public:
  /// Constructor.
  FluxLimiter(scalar* solution_vector, Hermes::vector<Space *> spaces, Hermes::vector<Solution *> solutions);

  /// Destructor.
   ~FluxLimiter();

  /// Do the limiting.
  void limit_according_to_detector(std::set<int>& discontinuous_elements);

protected:
  /// Members.
  Hermes::vector<Space *> spaces;
  Hermes::vector<Solution *> solutions;
  scalar* solution_vector;
};

FluxLimiter::FluxLimiter(scalar* solution_vector, Hermes::vector<Space *> spaces, Hermes::vector<Solution *> solutions) : solution_vector(solution_vector), spaces(spaces), 
  solutions(solutions)
{};

FluxLimiter::~FluxLimiter()
{};

void FluxLimiter::limit_according_to_detector(std::set<int>& discontinuous_elements)
{
  // First adjust the solution_vector.
  for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
      AsmList al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0)
          solution_vector[al.dof[shape_i]] = 0.0;
    }

  // Now adjust the solutions.
  Solution::vector_to_solutions(solution_vector, spaces, solutions);
};