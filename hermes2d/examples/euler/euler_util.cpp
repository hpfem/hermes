#include "euler_util.h"
#include "limits.h"

// Calculates energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa)
{
  return pressure/(kappa - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  return (kappa - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  return std::sqrt(kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy, kappa) / rho);
}

CFLCalculation::CFLCalculation(double CFL_number, double kappa) : CFL_number(CFL_number), kappa(kappa)
{
}

void CFLCalculation::calculate(Hermes::vector<Solution*> solutions, Mesh* mesh, double & time_step)
{
  // Create spaces of constant functions over the given mesh.
  L2Space constant_rho_space(mesh, 0);
  L2Space constant_rho_v_x_space(mesh, 0);
  L2Space constant_rho_v_y_space(mesh, 0);
  L2Space constant_energy_space(mesh, 0);

  scalar* sln_vector = new scalar[constant_rho_space.get_num_dofs() * 4];

  OGProjection::project_global(Hermes::vector<Space*>(&constant_rho_space, &constant_rho_v_x_space, &constant_rho_v_y_space, &constant_energy_space), solutions, sln_vector);

  // Determine the time step according to the CFL condition.

  double min_condition = 0;
  Element *e;
  for_all_active_elements(e, mesh) {
    AsmList al;
    constant_rho_space.get_element_assembly_list(e, &al);
    double rho = sln_vector[al.dof[0]];
    constant_rho_v_x_space.get_element_assembly_list(e, &al);
    double v1 = sln_vector[al.dof[0]] / rho;
    constant_rho_v_y_space.get_element_assembly_list(e, &al);
    double v2 = sln_vector[al.dof[0]] / rho;
    constant_energy_space.get_element_assembly_list(e, &al);
    double energy = sln_vector[al.dof[0]];
      
    double condition = e->get_area() * CFL_number / (std::sqrt(v1*v1 + v2*v2) + QuantityCalculator::calc_sound_speed(rho, rho*v1, rho*v2, energy, kappa));
      
    if(condition < min_condition || min_condition == 0.)
      min_condition = condition;
  }

  time_step = min_condition;

  delete [] sln_vector;
}

void CFLCalculation::calculate_semi_implicit(Hermes::vector<Solution*> solutions, Mesh* mesh, double & time_step)
{
  // Create spaces of constant functions over the given mesh.
  L2Space constant_rho_space(mesh, 0);
  L2Space constant_rho_v_x_space(mesh, 0);
  L2Space constant_rho_v_y_space(mesh, 0);
  L2Space constant_energy_space(mesh, 0);

  scalar* sln_vector = new scalar[constant_rho_space.get_num_dofs() * 4];

  OGProjection::project_global(Hermes::vector<Space*>(&constant_rho_space, &constant_rho_v_x_space, &constant_rho_v_y_space, &constant_energy_space), solutions, sln_vector);

  // Determine the time step according to the CFL condition.

  double min_condition = 0;
  Element *e;
  double w[4];
  for_all_active_elements(e, mesh) {
    AsmList al;
    constant_rho_space.get_element_assembly_list(e, &al);
    w[0] = sln_vector[al.dof[0]];
    constant_rho_v_x_space.get_element_assembly_list(e, &al);
    w[1] = sln_vector[al.dof[0]];
    constant_rho_v_y_space.get_element_assembly_list(e, &al);
    w[2] = sln_vector[al.dof[0]];
    constant_energy_space.get_element_assembly_list(e, &al);
    w[3] = sln_vector[al.dof[0]];
      
    double edge_length_max_lambda = 0.0;

    for(unsigned int edge_i = 0; edge_i < e->nvert; edge_i++) {
      // Initialization.      
      SurfPos surf_pos;
      surf_pos.marker = e->marker;
      surf_pos.surf_num = edge_i;
      int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 5);
      Geom<double>* geom = init_geom_surf(solutions[1]->get_refmap(), &surf_pos, eo);
      int np = solutions[1]->get_quad_2d()->get_num_points(eo);
      
      // Calculation of the edge length.
      double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));

      // Calculation of the maximum eigenvalue of the matrix P.
      double max_eigen_value = 0.0;
      for(int point_i = 0; point_i < np; point_i++) {
        // Transform to the local coordinates.
        double transformed[4];
        transformed[0] = w[0];
        transformed[1] = geom->nx[point_i] * w[1] + geom->ny[point_i] * w[2];
        transformed[2] = -geom->ny[point_i] * w[1] + geom->nx[point_i] * w[2];
        transformed[3] = w[3];

        // Calc sound speed.
        double a = QuantityCalculator::calc_sound_speed(transformed[0], transformed[1], transformed[2], transformed[3], kappa);
        
        // Calc max eigenvalue.
        if(transformed[1] / transformed[0] - a > max_eigen_value || point_i == 0)
          max_eigen_value = transformed[1] / transformed[0] - a;
        if(transformed[1] / transformed[0] > max_eigen_value)
          max_eigen_value = transformed[1] / transformed[0];
        if(transformed[1] / transformed[0] + a> max_eigen_value)
          max_eigen_value = transformed[1] / transformed[0] + a;
      }

      if(edge_length * max_eigen_value > edge_length_max_lambda || edge_i == 0)
        edge_length_max_lambda = edge_length * max_eigen_value;
      
      geom->free();
      delete geom;
    }

    
    double condition = e->get_area() * CFL_number / edge_length_max_lambda;
      
    if(condition < min_condition || min_condition == 0.)
      min_condition = condition;
  }

  time_step = min_condition;

  delete [] sln_vector;
}

void CFLCalculation::set_number(double new_CFL_number)
{
  this->CFL_number = new_CFL_number;
}

ADEStabilityCalculation::ADEStabilityCalculation(double AdvectionRelativeConstant, double DiffusionRelativeConstant, double epsilon) 
    : AdvectionRelativeConstant(AdvectionRelativeConstant), DiffusionRelativeConstant(DiffusionRelativeConstant), epsilon(epsilon)
{
}

void ADEStabilityCalculation::calculate(Hermes::vector<Solution*> solutions, Mesh* mesh, double & time_step)
{
  // Create spaces of constant functions over the given mesh.
  L2Space constant_rho_space(mesh, 0);
  L2Space constant_rho_v_x_space(mesh, 0);
  L2Space constant_rho_v_y_space(mesh, 0);

  scalar* sln_vector = new scalar[constant_rho_space.get_num_dofs() * 3];

  OGProjection::project_global(Hermes::vector<Space*>(&constant_rho_space, &constant_rho_v_x_space, &constant_rho_v_y_space), solutions, sln_vector);

  // Determine the time step according to the conditions.
  double min_condition_advection = 0.;
  double min_condition_diffusion = 0.;
  Element *e;
  for_all_active_elements(e, mesh) {
    AsmList al;
    constant_rho_space.get_element_assembly_list(e, &al);
    double rho = sln_vector[al.dof[0]];
    constant_rho_v_x_space.get_element_assembly_list(e, &al);
    double v1 = sln_vector[al.dof[0]] / rho;
    constant_rho_v_y_space.get_element_assembly_list(e, &al);
    double v2 = sln_vector[al.dof[0]] / rho;
      
    double condition_advection = AdvectionRelativeConstant * approximate_inscribed_circle_radius(e) / std::sqrt(v1*v1 + v2*v2);
    double condition_diffusion = DiffusionRelativeConstant * e->get_area() / epsilon;
      
    if(condition_advection < min_condition_advection || min_condition_advection == 0.)
      min_condition_advection = condition_advection;

    if(condition_diffusion < min_condition_diffusion || min_condition_diffusion == 0.)
      min_condition_diffusion = condition_diffusion;
  }

  time_step = std::min(min_condition_advection, min_condition_diffusion);

  delete [] sln_vector;
}

double ADEStabilityCalculation::approximate_inscribed_circle_radius(Element * e)
{
  double h = std::sqrt(std::pow(e->vn[(0 + 1) % e->get_num_surf()]->x - e->vn[0]->x, 2) + std::pow(e->vn[(0 + 1) % e->get_num_surf()]->y - e->vn[0]->y, 2));
    for(int edge_i = 0; edge_i < e->get_num_surf(); edge_i++) {
      double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));
      if(edge_length < h)
        h = edge_length;
    }
  return h / 2;
}

DiscontinuityDetector::DiscontinuityDetector(Hermes::vector<Space *> spaces, 
                        Hermes::vector<Solution *> solutions) : spaces(spaces), solutions(solutions)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      error("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

DiscontinuityDetector::~DiscontinuityDetector()
{};

double DiscontinuityDetector::calculate_h(Element* e, int polynomial_order)
{
  double h = std::sqrt(std::pow(e->vn[(0 + 1) % e->get_num_surf()]->x - e->vn[0]->x, 2) + std::pow(e->vn[(0 + 1) % e->get_num_surf()]->y - e->vn[0]->y, 2));
  for(int edge_i = 0; edge_i < e->get_num_surf(); edge_i++) {
    double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));
    if(edge_length < h)
      h = edge_length;
  }
  return std::pow(h, (0.5 * (H2D_GET_H_ORDER(spaces[0]->get_element_order(e->id)) 
                                                     + 
                                                     H2D_GET_V_ORDER(spaces[0]->get_element_order(e->id)))
                                              + 1) / 2);
}

std::set<int>& DiscontinuityDetector::get_discontinuous_element_ids(double threshold)
{
  Element* e;
  for_all_active_elements(e, mesh) {
    bool element_inserted = false;
    for(int edge_i = 0; edge_i < e->get_num_surf() && !element_inserted; edge_i++)
      if(calculate_relative_flow_direction(e, edge_i) < 0 && !e->en[edge_i]->bnd) {
        double jumps[4];
        calculate_jumps(e, edge_i, jumps);
        double diameter_indicator = calculate_h(e, spaces[0]->get_element_order(e->id));
        double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_num_surf()]->y - e->vn[edge_i]->y, 2));
        double norms[4];
        calculate_norms(e, edge_i, norms);

        // Number of component jumps tested.
        unsigned int component_checked_number = 4;

        for(unsigned int component_i = 0; component_i < component_checked_number; component_i++) {
          if(norms[component_i] < 1E-8)
            continue;
          double discontinuity_detector = jumps[component_i] / (diameter_indicator * edge_length * norms[component_i]);
          if(discontinuity_detector > threshold) {
            discontinuous_element_ids.insert(e->id);
            element_inserted = true;
            break;
          }
        }
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

  geom->free();
  delete geom;
  delete [] jwt;
  density_vel_x->free_fn();
  density_vel_y->free_fn();
  delete density_vel_x;
  delete density_vel_y;

  return result;
};

void DiscontinuityDetector::calculate_jumps(Element* e, int edge_i, double result[4])
{
  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8);
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo);
  int np = solutions[0]->get_quad_2d()->get_num_points(eo);

  // Initialize the NeighborSearch.
  NeighborSearch ns(e, mesh);
  ns.set_active_edge(edge_i);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Go through all neighbors.
  for(int neighbor_i = 0; neighbor_i < ns.get_num_neighbors(); neighbor_i++) {
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
    Func<scalar>* energy = init_fn(solutions[3], eo);

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
    Func<scalar>* energy_neighbor = init_fn(solutions[3], eo);

    DiscontinuousFunc<scalar> density_discontinuous(density, density_neighbor, true);
    DiscontinuousFunc<scalar> density_vel_x_discontinuous(density_vel_x, density_vel_x_neighbor, true);
    DiscontinuousFunc<scalar> density_vel_y_discontinuous(density_vel_y, density_vel_y_neighbor, true);
    DiscontinuousFunc<scalar> energy_discontinuous(energy, energy_neighbor, true);

    for(int point_i = 0; point_i < np; point_i++) {
      result[0] += jwt[point_i] * (density_discontinuous.get_val_central(point_i) - density_discontinuous.get_val_neighbor(point_i)); 
      result[1] += jwt[point_i] * (density_vel_x_discontinuous.get_val_central(point_i) - density_vel_x_discontinuous.get_val_neighbor(point_i));
      result[2] += jwt[point_i] * (density_vel_y_discontinuous.get_val_central(point_i) - density_vel_y_discontinuous.get_val_neighbor(point_i));
      result[3] += jwt[point_i] * (energy_discontinuous.get_val_central(point_i) - energy_discontinuous.get_val_neighbor(point_i));
    }
    
    geom->free();
    delete geom;
    delete [] jwt;
    density->free_fn();
    density_vel_x->free_fn();
    density_vel_y->free_fn();
    energy->free_fn();
    density_neighbor->free_fn();
    density_vel_x_neighbor->free_fn();
    density_vel_y_neighbor->free_fn();
    energy_neighbor->free_fn();
    
    delete density;
    delete density_vel_x;
    delete density_vel_y;
    delete energy;
    delete density_neighbor;
    delete density_vel_x_neighbor;
    delete density_vel_y_neighbor;
    delete energy_neighbor;
  }

  result[0] = std::abs(result[0]);
  result[1] = std::abs(result[1]);
  result[2] = std::abs(result[2]);
  result[3] = std::abs(result[3]);
};

void DiscontinuityDetector::calculate_norms(Element* e, int edge_i, double result[4])
{
  // Set active element to the solutions.
  solutions[0]->set_active_element(e);
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);
  solutions[3]->set_active_element(e);
  
  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8);
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
  Func<scalar>* energy = init_fn(solutions[3], eo);

  for(int point_i = 0; point_i < np; point_i++) {
    result[0] = std::max(result[0], std::abs(density->val[point_i]));
    result[1] = std::max(result[1], std::abs(density_vel_x->val[point_i]));
    result[2] = std::max(result[2], std::abs(density_vel_y->val[point_i]));
    result[3] = std::max(result[3], std::abs(energy->val[point_i]));
  }

  geom->free();
  delete geom;
  delete [] jwt;
  
  density->free_fn();
  density_vel_x->free_fn();
  density_vel_y->free_fn();
  energy->free_fn();
    
  delete density;
  delete density_vel_x;
  delete density_vel_y;
  delete energy;
};

FluxLimiter::FluxLimiter(scalar* solution_vector, Hermes::vector<Space *> spaces, Hermes::vector<Solution *> solutions) : solution_vector(solution_vector), spaces(spaces), 
  solutions(solutions)
{};

FluxLimiter::~FluxLimiter()
{};

void FluxLimiter::limit_according_to_detector(std::set<int>& discontinuous_elements, Hermes::vector<Space *> coarse_spaces)
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

  if(coarse_spaces != Hermes::vector<Space *>()) {
    // Now set the element order to zero.
    Element* e;

    for_all_elements(e, spaces[0]->get_mesh())
      e->visited = false;

    for(unsigned int space_i = 0; space_i < spaces.size(); space_i++) {
      for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
        AsmList al;
        spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
        for(unsigned int shape_i = 0; shape_i < al.cnt; shape_i++) {
          if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.idx[shape_i])) > 0) {
            spaces[space_i]->get_mesh()->get_element(*it)->visited = true;
            bool all_sons_visited = true;
            for(unsigned int son_i = 0; son_i < 4; son_i++)
              if(!spaces[space_i]->get_mesh()->get_element(*it)->parent->sons[son_i]->visited) {
                all_sons_visited = false;
                break;
              }
            if(all_sons_visited)
              coarse_spaces[space_i]->set_element_order_internal(spaces[space_i]->get_mesh()->get_element(*it)->parent->id, 0);
            break;
          }
        }
      }
    }

    Space::assign_dofs(coarse_spaces);
  }
};

void MachNumberFilter::filter_fn(int n, Hermes::vector<scalar*> values, scalar* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(kappa * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / values.at(0)[i]);
}

void PressureFilter::filter_fn(int n, Hermes::vector<scalar*> values, scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (kappa - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}


void EntropyFilter::filter_fn(int n, Hermes::vector<scalar*> values, scalar* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / p_ext)
      / pow((values.at(0)[i] / rho_ext), kappa));
}