#include "euler_util.h"
#include "limits.h"
#include <limits>

// Calculates energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa)
{
  double to_return = pressure/(kappa - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = (kappa - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = std::sqrt(kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy, kappa) / rho);
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

CFLCalculation::CFLCalculation(double CFL_number, double kappa) : CFL_number(CFL_number), kappa(kappa)
{
}

void CFLCalculation::calculate(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step) const
{
  // Create spaces of constant functions over the given mesh->
  SpaceSharedPtr<double> constant_rho_space(new L2Space<double> (mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_x_space(new L2Space<double> (mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_y_space(new L2Space<double> (mesh, 0));
  SpaceSharedPtr<double> constant_energy_space(new L2Space<double> (mesh, 0));

  double* sln_vector = new double[constant_rho_space->get_num_dofs() * 4];

  OGProjection<double> ogProjection;
  ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double>  >(constant_rho_space, constant_rho_v_x_space, constant_rho_v_y_space, constant_energy_space), solutions, sln_vector);

  // Determine the time step according to the CFL condition.

  double min_condition = 0;
  Element *e;
  for_all_active_elements(e, mesh)
  {
    AsmList<double> al;
    constant_rho_space->get_element_assembly_list(e, &al);
    double rho = sln_vector[al.get_dof()[0]];
    constant_rho_v_x_space->get_element_assembly_list(e, &al);
    double v1 = sln_vector[al.get_dof()[0]] / rho;
    constant_rho_v_y_space->get_element_assembly_list(e, &al);
    double v2 = sln_vector[al.get_dof()[0]] / rho;
    constant_energy_space->get_element_assembly_list(e, &al);
    double energy = sln_vector[al.get_dof()[0]];

    e->calc_area();
    double condition = e->area * CFL_number / (std::sqrt(v1*v1 + v2*v2) + QuantityCalculator::calc_sound_speed(rho, rho*v1, rho*v2, energy, kappa));

    if(condition < min_condition || min_condition == 0.)
      min_condition = condition;
  }

  time_step = min_condition;

  delete [] sln_vector;
}

void CFLCalculation::calculate_semi_implicit(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step) const
{
  // Create spaces of constant functions over the given mesh->
  SpaceSharedPtr<double> constant_rho_space(new L2Space<double>(mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_x_space(new L2Space<double>(mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_y_space(new L2Space<double>(mesh, 0));
  SpaceSharedPtr<double> constant_energy_space(new L2Space<double>(mesh, 0));

  double* sln_vector = new double[constant_rho_space->get_num_dofs() * 4];

  OGProjection<double> ogProjection;
  ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double>  >(constant_rho_space, constant_rho_v_x_space, constant_rho_v_y_space, constant_energy_space), solutions, sln_vector);

  // Determine the time step according to the CFL condition.

  double min_condition = 0;
  Element *e;
  double w[4];
  for_all_active_elements(e, mesh)
  {
    AsmList<double> al;
    constant_rho_space->get_element_assembly_list(e, &al);
    w[0] = sln_vector[al.get_dof()[0]];
    constant_rho_v_x_space->get_element_assembly_list(e, &al);
    w[1] = sln_vector[al.get_dof()[0]];
    constant_rho_v_y_space->get_element_assembly_list(e, &al);
    w[2] = sln_vector[al.get_dof()[0]];
    constant_energy_space->get_element_assembly_list(e, &al);
    w[3] = sln_vector[al.get_dof()[0]];

    double edge_length_max_lambda = 0.0;

    solutions[0]->set_active_element(e);
    for(unsigned int edge_i = 0; edge_i < e->get_nvert(); edge_i++) {
      // Initialization.
      SurfPos surf_pos;
      surf_pos.marker = e->marker;
      surf_pos.surf_num = edge_i;
      int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 20, e->get_mode());
      double3* tan = NULL;
      Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);
      int np = solutions[1]->get_quad_2d()->get_num_points(eo, e->get_mode());

      // Calculation of the edge length.
      double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->y - e->vn[edge_i]->y, 2));

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

    e->calc_area();
    double condition = e->area * CFL_number / edge_length_max_lambda;

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

void ADEStabilityCalculation::calculate(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step)
{
  // Create spaces of constant functions over the given mesh->
  SpaceSharedPtr<double> constant_rho_space(new L2Space<double>(mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_x_space(new L2Space<double>(mesh, 0));
  SpaceSharedPtr<double> constant_rho_v_y_space(new L2Space<double>(mesh, 0));

  double* sln_vector = new double[constant_rho_space->get_num_dofs() * 3];

  OGProjection<double> ogProjection;
  ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double>  >(constant_rho_space, constant_rho_v_x_space, constant_rho_v_y_space), solutions, sln_vector);

  // Determine the time step according to the conditions.
  double min_condition_advection = 0.;
  double min_condition_diffusion = 0.;
  Element *e;
  for_all_active_elements(e, mesh)
  {
    AsmList<double> al;
    constant_rho_space->get_element_assembly_list(e, &al);
    double rho = sln_vector[al.get_dof()[0]];
    constant_rho_v_x_space->get_element_assembly_list(e, &al);
    double v1 = sln_vector[al.get_dof()[0] + constant_rho_space->get_num_dofs()] / rho;
    constant_rho_v_y_space->get_element_assembly_list(e, &al);
    double v2 = sln_vector[al.get_dof()[0] + 2 * constant_rho_space->get_num_dofs()] / rho;

    e->calc_area();
    e->calc_diameter();
    double condition_advection = AdvectionRelativeConstant * e->diameter / std::sqrt(v1*v1 + v2*v2);
    double condition_diffusion = DiffusionRelativeConstant * e->area / epsilon;

    if(condition_advection < min_condition_advection || min_condition_advection == 0.)
      min_condition_advection = condition_advection;

    if(condition_diffusion < min_condition_diffusion || min_condition_diffusion == 0.)
      min_condition_diffusion = condition_diffusion;
  }

  time_step = std::min(min_condition_advection, min_condition_diffusion);

  delete [] sln_vector;
}

DiscontinuityDetector::DiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double>  > spaces, 
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : spaces(spaces), solutions(solutions)
{
  for(int i = 0; i < solutions.size(); i++)
    this->solutionsInternal.push_back((Solution<double>*)solutions[i].get());
};

DiscontinuityDetector::~DiscontinuityDetector()
{};

KrivodonovaDiscontinuityDetector::KrivodonovaDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double>  > spaces, 
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : DiscontinuityDetector(spaces, solutions)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      throw Hermes::Exceptions::Exception("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

KrivodonovaDiscontinuityDetector::~KrivodonovaDiscontinuityDetector()
{};

double KrivodonovaDiscontinuityDetector::calculate_h(Element* e, int polynomial_order)
{
  double h = std::sqrt(std::pow(e->vn[(0 + 1) % e->get_nvert()]->x - e->vn[0]->x, 2) + std::pow(e->vn[(0 + 1) % e->get_nvert()]->y - e->vn[0]->y, 2));
  for(int edge_i = 0; edge_i < e->get_nvert(); edge_i++) {
    double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->y - e->vn[edge_i]->y, 2));
    if(edge_length < h)
      h = edge_length;
  }
  return std::pow(h, (0.5 * (H2D_GET_H_ORDER(spaces[0]->get_element_order(e->id)) 
    + 
    H2D_GET_V_ORDER(spaces[0]->get_element_order(e->id)))
    + 1) / 2);
}

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids()
{
  return get_discontinuous_element_ids(1.0);
};

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids(double threshold)
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    bool element_inserted = false;
    for(int edge_i = 0; edge_i < e->get_nvert() && !element_inserted; edge_i++)
      if(calculate_relative_flow_direction(e, edge_i) < 0 && !e->en[edge_i]->bnd)
      {
        double jumps[4];
        calculate_jumps(e, edge_i, jumps);
        double diameter_indicator = calculate_h(e, spaces[0]->get_element_order(e->id));
        double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->y - e->vn[edge_i]->y, 2));
        double norms[4];
        calculate_norms(e, edge_i, norms);

        // Number of component jumps tested.
        unsigned int component_checked_number = 1;
        for(unsigned int component_i = 0; component_i < component_checked_number; component_i++) {
          if(norms[component_i] < 1E-8)
            continue;
          double discontinuity_detector = jumps[component_i] / (diameter_indicator * edge_length * norms[component_i]);
          if(discontinuity_detector > threshold)
          {
            discontinuous_element_ids.insert(e->id);
            element_inserted = true;
            break;
          }
        }
      }
  }
  return discontinuous_element_ids;
};

double KrivodonovaDiscontinuityDetector::calculate_relative_flow_direction(Element* e, int edge_i)
{
  // Set active element to the two solutions (density_vel_x, density_vel_y).
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 20, e->get_mode());
  double3* pt = solutions[1]->get_quad_2d()->get_points(eo, e->get_mode());
  int np = solutions[1]->get_quad_2d()->get_num_points(eo, e->get_mode());

  double3* tan;
  Geom<double>* geom = init_geom_surf(solutions[1]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);

  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Calculate.
  Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
  Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);

  double result = 0.0;
  for(int point_i = 0; point_i < np; point_i++)
    result += jwt[point_i] * density_vel_x->val[point_i] * geom->nx[point_i] + density_vel_y->val[point_i] * geom->ny[point_i];

  geom->free();
  delete geom;
  delete [] jwt;
  delete density_vel_x;
  delete density_vel_y;

  return result;
};

void KrivodonovaDiscontinuityDetector::calculate_jumps(Element* e, int edge_i, double result[4])
{
  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8, e->get_mode());
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo, e->get_mode());
  int np = solutions[0]->get_quad_2d()->get_num_points(eo, e->get_mode());

  // Initialize the NeighborSearch.
  NeighborSearch<double> ns(e, mesh);
  ns.set_active_edge(edge_i);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Go through all neighbors.
  for(int neighbor_i = 0; neighbor_i < ns.get_num_neighbors(); neighbor_i++) {
    ns.set_active_segment(neighbor_i);

    // Set active element to the solutions.
    solutions[0]->set_active_element(e);
    solutions[1]->set_active_element(e);
    solutions[2]->set_active_element(e);
    solutions[3]->set_active_element(e);

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_central_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
    }

    double3* tan;
    Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);
    double* jwt = new double[np];
    for(int i = 0; i < np; i++)
      jwt[i] = pt[i][2] * tan[i][2];

    // Prepare functions on the central element.
    Func<double>* density = init_fn(solutions[0].get(), eo);
    Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
    Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);
    Func<double>* energy = init_fn(solutions[3].get(), eo);

    // Set neighbor element to the solutions.
    solutions[0]->set_active_element(ns.get_neighb_el());
    solutions[1]->set_active_element(ns.get_neighb_el());
    solutions[2]->set_active_element(ns.get_neighb_el());
    solutions[3]->set_active_element(ns.get_neighb_el());

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_neighbor_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
    }

    // Prepare functions on the neighbor element.
    Func<double>* density_neighbor = init_fn(solutions[0].get(), eo);
    Func<double>* density_vel_x_neighbor = init_fn(solutions[1].get(), eo);
    Func<double>* density_vel_y_neighbor = init_fn(solutions[2].get(), eo);
    Func<double>* energy_neighbor = init_fn(solutions[3].get(), eo);

    DiscontinuousFunc<double> density_discontinuous(density, density_neighbor, true);
    DiscontinuousFunc<double> density_vel_x_discontinuous(density_vel_x, density_vel_x_neighbor, true);
    DiscontinuousFunc<double> density_vel_y_discontinuous(density_vel_y, density_vel_y_neighbor, true);
    DiscontinuousFunc<double> energy_discontinuous(energy, energy_neighbor, true);

    for(int point_i = 0; point_i < np; point_i++) {
      result[0] += jwt[point_i] * std::abs(density_discontinuous.val[point_i] - density_discontinuous.val_neighbor[point_i]); 
      result[1] += jwt[point_i] * std::abs(density_vel_x_discontinuous.val[point_i] - density_vel_x_discontinuous.val_neighbor[point_i]);
      result[2] += jwt[point_i] * std::abs(density_vel_y_discontinuous.val[point_i] - density_vel_y_discontinuous.val_neighbor[point_i]);
      result[3] += jwt[point_i] * std::abs(energy_discontinuous.val[point_i] - energy_discontinuous.val_neighbor[point_i]);
    }

    geom->free();
    delete geom;
    delete [] jwt;

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

void KrivodonovaDiscontinuityDetector::calculate_norms(Element* e, int edge_i, double result[4])
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

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8, e->get_mode());
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo, e->get_mode());
  int np = solutions[0]->get_quad_2d()->get_num_points(eo, e->get_mode());

  double3* tan;
  Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Calculate.
  Func<double>* density = init_fn(solutions[0].get(), eo);
  Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
  Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);
  Func<double>* energy = init_fn(solutions[3].get(), eo);

  for(int point_i = 0; point_i < np; point_i++) {
    result[0] = std::max(result[0], std::abs(density->val[point_i]));
    result[1] = std::max(result[1], std::abs(density_vel_x->val[point_i]));
    result[2] = std::max(result[2], std::abs(density_vel_y->val[point_i]));
    result[3] = std::max(result[3], std::abs(energy->val[point_i]));
  }

  geom->free();
  delete geom;
  delete [] jwt;

  delete density;
  delete density_vel_x;
  delete density_vel_y;
  delete energy;
};

KuzminDiscontinuityDetector::KuzminDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double>  > spaces, 
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions, bool limit_all_orders_independently) : DiscontinuityDetector(spaces, solutions), limit_all_orders_independently(limit_all_orders_independently)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      throw Hermes::Exceptions::Exception("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

KuzminDiscontinuityDetector::~KuzminDiscontinuityDetector()
{};

std::set<int>& KuzminDiscontinuityDetector::get_discontinuous_element_ids()
{
  Element* e;
  discontinuous_element_ids.clear();
  oscillatory_element_idsRho.clear();
  oscillatory_element_idsRhoVX.clear();
  oscillatory_element_idsRhoVY.clear();
  oscillatory_element_idsRhoE.clear();
  for_all_active_elements(e, mesh)
  {
    if(!limit_all_orders_independently)
      if(this->second_order_discontinuous_element_ids.find(e->id) == this->second_order_discontinuous_element_ids.end())
        continue;
    if(e->get_nvert() == 3)
      throw Hermes::Exceptions::Exception("So far this limiter is implemented just for quads.");
    double u_c[4], u_dx_c[4], u_dy_c[4];

    double x_center, y_center, x_center_ref, y_center_ref;
    e->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(e, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_values(e, u_c, x_center_ref, y_center_ref);

    // Vertex values.
    double u_i[4][4];
    find_vertex_values(e, u_i);

    // Boundaries for alpha_i calculation.
    double u_i_min_first_order[4][4];
    double u_i_max_first_order[4][4];
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
      {
        u_i_min_first_order[i][j] = std::numeric_limits<double>::infinity();
        u_i_max_first_order[i][j] = -std::numeric_limits<double>::infinity();
      }
      find_u_i_min_max_first_order(e, u_i_min_first_order, u_i_max_first_order);

      // alpha_i calculation.
      double alpha_i_first_order[4];
      find_alpha_i_first_order(u_i_min_first_order, u_i_max_first_order, u_c, u_i, alpha_i_first_order);

      // measure.
      for(unsigned int i = 0; i < 4; i++)
      {
        if(1.0 > alpha_i_first_order[i])
        {
          // check for sanity.
          if(std::abs(u_c[i]) > 1E-12)
            discontinuous_element_ids.insert(e->id);
        }
        if(std::abs(u_c[i]) < 1e-3)
          continue;
        bool bnd = false;
        for(unsigned int j = 0; j < 4; j++)
          if(e->en[j]->bnd)
            bnd = true;
        if(bnd)
          continue;
        double high_limit = -std::numeric_limits<double>::infinity();
        double low_limit = std::numeric_limits<double>::infinity();
        for(unsigned int j = 0; j < 4; j++)
        {
          if(u_i_max_first_order[i][j] > high_limit)
            high_limit = u_i_max_first_order[i][j];
          if(u_i_min_first_order[i][j] < low_limit)
            low_limit = u_i_min_first_order[i][j];
        }
        if(u_c[i] > high_limit && std::abs(u_c[i] - high_limit) > 1e-3)
        {
          high_limit = (u_i_max_first_order[i][0] + u_i_max_first_order[i][1] + u_i_max_first_order[i][2] + u_i_max_first_order[i][3] + u_i_min_first_order[i][0] + u_i_min_first_order[i][1] + u_i_min_first_order[i][2] + u_i_min_first_order[i][3]) / 8.0;
          switch(i)
          {
          case 0:
            oscillatory_element_idsRho.insert(std::pair<int, double>(e->id, high_limit));
            break;
          case 1:
            oscillatory_element_idsRhoVX.insert(std::pair<int, double>(e->id, high_limit));
            break;
          case 2:
            oscillatory_element_idsRhoVY.insert(std::pair<int, double>(e->id, high_limit));
            break;
          case 3:
            oscillatory_element_idsRhoE.insert(std::pair<int, double>(e->id, high_limit));
            break;
          }
        }
        if(u_c[i] < low_limit && std::abs(u_c[i] - low_limit) > 1e-3)
        {
          low_limit = (u_i_max_first_order[i][0] + u_i_max_first_order[i][1] + u_i_max_first_order[i][2] + u_i_max_first_order[i][3] + u_i_min_first_order[i][0] + u_i_min_first_order[i][1] + u_i_min_first_order[i][2] + u_i_min_first_order[i][3]) / 8.0;
          switch(i)
          {
          case 0:
            oscillatory_element_idsRho.insert(std::pair<int, double>(e->id, low_limit));
            break;
          case 1:
            oscillatory_element_idsRhoVX.insert(std::pair<int, double>(e->id, low_limit));
            break;
          case 2:
            oscillatory_element_idsRhoVY.insert(std::pair<int, double>(e->id, low_limit));
            break;
          case 3:
            oscillatory_element_idsRhoE.insert(std::pair<int, double>(e->id, low_limit));
            break;
          }
        }
      }
  }

  return discontinuous_element_ids;
}

std::set<int>& KuzminDiscontinuityDetector::get_second_order_discontinuous_element_ids()
{
  Element* e;

  for_all_active_elements(e, mesh)
  {
    if(e->get_nvert() == 3)
      throw Hermes::Exceptions::Exception("So far this limiter is implemented just for quads.");
    double c_x, c_y;
    e->get_center(c_x, c_y);

    double*** values = new double**[this->solutions.size()];
    for(int i = 0; i < this->solutions.size(); i++)
      values[i] = solutionsInternal[i]->get_ref_values_transformed(e, c_x, c_y);

    // Vertex values.
    double u_d_i[4][4][2];
    find_vertex_derivatives(e, u_d_i);

    // Boundaries for alpha_i calculation.
    double u_d_i_min_second_order[4][4][2];
    double u_d_i_max_second_order[4][4][2];
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        for(int k = 0; k < 2; k++)
        {
          u_d_i_min_second_order[i][j][k] = std::numeric_limits<double>::infinity();
          u_d_i_max_second_order[i][j][k] = -std::numeric_limits<double>::infinity();
        }

        find_u_i_min_max_second_order(e, u_d_i_min_second_order, u_d_i_max_second_order);

        // alpha_i calculation.
        double alpha_i_second_order[4];
        find_alpha_i_second_order(u_d_i_min_second_order, u_d_i_max_second_order, values, values, u_d_i, alpha_i_second_order);

        // measure.
        for(unsigned int i = 0; i < 4; i++)
          if(1.0 > alpha_i_second_order[i])
          {
            // check for sanity.
            if(std::abs(values[i][0][1]) > 1E-12 || std::abs(values[i][0][2]) > 1E-12)
              second_order_discontinuous_element_ids.insert(e->id);
          }
  }

  return second_order_discontinuous_element_ids;
}

bool KuzminDiscontinuityDetector::get_limit_all_orders_independently()
{
  return this->limit_all_orders_independently;
}

void KuzminDiscontinuityDetector::find_centroid_values(Hermes::Hermes2D::Element* e, double u_c[4], double x_ref, double y_ref)
{
  for(unsigned int i = 0; i < this->solutions.size(); i++)
    u_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, x_ref, y_ref, 0, 0);
}

void KuzminDiscontinuityDetector::find_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dx_c[4], double u_dy_c[4], double x_ref, double y_ref)
{
  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    u_dx_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, x_ref, y_ref, 0, 1);
    u_dy_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, x_ref, y_ref, 0, 2);
  }
}

void KuzminDiscontinuityDetector::find_second_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4])
{
  double c_x, c_y;
  double c_ref_x, c_ref_y;
  if(e->get_nvert() == 3)
  {
    c_x = (0.33333333333333333) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
    c_y = (0.33333333333333333) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
    c_x = (0.25) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
    c_y = (0.25) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    solutions[i]->set_active_element(e);
    solutions[i]->get_refmap()->untransform(e, c_x, c_y, c_ref_x, c_ref_y);
    u_dxx_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 3);
    u_dyy_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 4);
    u_dxy_c[i] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 5);
  }
}

void KuzminDiscontinuityDetector::find_vertex_values(Hermes::Hermes2D::Element* e, double vertex_values[4][4])
{
  double c_ref_x, c_ref_y;
  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    for(unsigned int j = 0; j < e->get_nvert(); j++)
    {
      solutions[i]->get_refmap()->set_active_element(e);
      solutions[i]->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_values[i][j] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 0);
    }
  }
}

void KuzminDiscontinuityDetector::find_vertex_derivatives(Hermes::Hermes2D::Element* e, double vertex_derivatives[4][4][2])
{
  double c_ref_x, c_ref_y;
  for(unsigned int i = 0; i < this->solutions.size(); i++)
  {
    for(unsigned int j = 0; j < e->get_nvert(); j++)
    {
      solutions[i]->get_refmap()->set_active_element(e);
      solutions[i]->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_derivatives[i][j][0] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 1);
      vertex_derivatives[i][j][1] = solutionsInternal[i]->get_ref_value_transformed(e, c_ref_x, c_ref_y, 0, 2);
    }
  }
}

void KuzminDiscontinuityDetector::find_u_i_min_max_first_order(Hermes::Hermes2D::Element* e, double u_i_min[4][4], double u_i_max[4][4])
{
  for(unsigned int j = 0; j < e->get_nvert(); j++)
  {
    Hermes::Hermes2D::NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;
    ns.set_active_edge(j);

    // First beginning neighbors on every edge.
    double u_c[4];
    ns.set_active_segment(0);

    Element* en = ns.get_neighb_el();
    double x_center, y_center, x_center_ref, y_center_ref;
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_values(en, u_c, x_center_ref, y_center_ref);

    for(unsigned int min_i = 0; min_i < 4; min_i++)
      if(u_i_min[min_i][j] > u_c[min_i])
        u_i_min[min_i][j] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 4; max_i++)
      if(u_i_max[max_i][j] < u_c[max_i])
        u_i_max[max_i][j] = u_c[max_i];

    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);

    en = ns.get_neighb_el();
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_values(en, u_c, x_center_ref, y_center_ref);

    for(unsigned int min_i = 0; min_i < 4; min_i++)
      if(u_i_min[min_i][(j + 1) % e->get_nvert()] > u_c[min_i])
        u_i_min[min_i][(j + 1) % e->get_nvert()] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 4; max_i++)
      if(u_i_max[max_i][(j + 1) % e->get_nvert()] < u_c[max_i])
        u_i_max[max_i][(j + 1) % e->get_nvert()] = u_c[max_i];

    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    Hermes::Hermes2D::NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert());
    ns_1.set_active_segment(0);

    en = ns_1.get_neighb_el();
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_values(en, u_c, x_center_ref, y_center_ref);

    for(unsigned int min_i = 0; min_i < 4; min_i++)
      if(u_i_min[min_i][j] > u_c[min_i])
        u_i_min[min_i][j] = u_c[min_i];
    for(unsigned int max_i = 0; max_i < 4; max_i++)
      if(u_i_max[max_i][j] < u_c[max_i])
        u_i_max[max_i][j] = u_c[max_i];
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_first_order(double u_i_min[4][4], double u_i_max[4][4], double u_c[4], double u_i[4][4], double alpha_i[4])
{
  for(unsigned int sol_i = 0; sol_i < 4; sol_i++)
  {
    alpha_i[sol_i] = 1;
    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
      // Sanity checks.
      if(std::abs(u_i[sol_i][vertex_i] - u_c[sol_i]) < 1E-6)
        continue;
      if(std::abs((u_i_min[sol_i][vertex_i] - u_c[sol_i]) / u_c[sol_i]) > 10)
        continue;
      if(std::abs((u_i_max[sol_i][vertex_i] - u_c[sol_i]) / u_c[sol_i]) > 10)
        continue;

      if(u_i[sol_i][vertex_i] < u_c[sol_i])
      {
        if((u_i_min[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_i_min[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]);
      }
      else
      {
        if((u_i_max[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_i_max[sol_i][vertex_i] - u_c[sol_i]) / (u_i[sol_i][vertex_i] - u_c[sol_i]);
      }
    }
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_first_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4], double u_c[4], double u_dx_c[4], double u_dy_c[4], double alpha_i_real[4])
{
  double c_x, c_y;
  if(e->get_nvert() == 3)
  {
    c_x = (1/3) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
    c_y = (1/3) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
    c_x = (1/4) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
    c_y = (1/4) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  alpha_i_real[0] = alpha_i_real[1] = alpha_i_real[2] = alpha_i_real[3] = -std::numeric_limits<double>::infinity();
  for(unsigned int sol_i = 0; sol_i < this->solutions.size(); sol_i++)
    for(unsigned int vertex_i = 0; vertex_i < e->get_nvert(); vertex_i++)
      if( (u_i[sol_i][vertex_i] - u_c[sol_i]) / (u_dx_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dy_c[sol_i] * (e->vn[vertex_i]->y - c_y)) > alpha_i_real[sol_i])
        alpha_i_real[sol_i] = (u_i[sol_i][vertex_i] - u_c[sol_i]) / (u_dx_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dy_c[sol_i] * (e->vn[vertex_i]->y - c_y));
}

void KuzminDiscontinuityDetector::find_u_i_min_max_second_order(Hermes::Hermes2D::Element* e, double u_d_i_min[4][4][2], double u_d_i_max[4][4][2])
{
  for(unsigned int j = 0; j < e->get_nvert(); j++)
  {
    Hermes::Hermes2D::NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;

    ns.set_active_edge(j);

    // First beginning neighbors on every edge.
    double u_dx_c[4], u_dy_c[4];
    ns.set_active_segment(0);

    Element* en = ns.get_neighb_el();
    double x_center, y_center, x_center_ref, y_center_ref;
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_derivatives(en, u_dx_c, u_dy_c, x_center_ref, y_center_ref);
    for(unsigned int min_i = 0; min_i < 4; min_i++)
    {
      if(u_d_i_min[min_i][j][0] > u_dx_c[min_i])
        u_d_i_min[min_i][j][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][j][1] > u_dy_c[min_i])
        u_d_i_min[min_i][j][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 4; max_i++)
    {
      if(u_d_i_max[max_i][j][0] < u_dx_c[max_i])
        u_d_i_max[max_i][j][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][j][1] < u_dy_c[max_i])
        u_d_i_max[max_i][j][1] = u_dy_c[max_i];
    }
    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);

    en = ns.get_neighb_el();
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_derivatives(en, u_dx_c, u_dy_c, x_center_ref, y_center_ref);

    for(unsigned int min_i = 0; min_i < 4; min_i++)
    {
      if(u_d_i_min[min_i][(j + 1) % e->get_nvert()][0] > u_dx_c[min_i])
        u_d_i_min[min_i][(j + 1) % e->get_nvert()][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][(j + 1) % e->get_nvert()][1] > u_dy_c[min_i])
        u_d_i_min[min_i][(j + 1) % e->get_nvert()][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 4; max_i++)
    {
      if(u_d_i_max[max_i][(j + 1) % e->get_nvert()][0] < u_dx_c[max_i])
        u_d_i_max[max_i][(j + 1) % e->get_nvert()][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][(j + 1) % e->get_nvert()][1] < u_dy_c[max_i])
        u_d_i_max[max_i][(j + 1) % e->get_nvert()][1] = u_dy_c[max_i];
    }

    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    Hermes::Hermes2D::NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert());
    ns_1.set_active_segment(0);

    en = ns_1.get_neighb_el();
    en->get_center(x_center, y_center);
    solutions[0]->get_refmap()->untransform(en, x_center, y_center, x_center_ref, y_center_ref);

    find_centroid_derivatives(en, u_dx_c, u_dy_c, x_center_ref, y_center_ref);

    for(unsigned int min_i = 0; min_i < 4; min_i++)
    {
      if(u_d_i_min[min_i][j][0] > u_dx_c[min_i])
        u_d_i_min[min_i][j][0] = u_dx_c[min_i];
      if(u_d_i_min[min_i][j][1] > u_dy_c[min_i])
        u_d_i_min[min_i][j][1] = u_dy_c[min_i];
    }
    for(unsigned int max_i = 0; max_i < 4; max_i++)
    {
      if(u_d_i_max[max_i][j][0] < u_dx_c[max_i])
        u_d_i_max[max_i][j][0] = u_dx_c[max_i];
      if(u_d_i_max[max_i][j][1] < u_dy_c[max_i])
        u_d_i_max[max_i][j][1] = u_dy_c[max_i];
    }
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_second_order(double u_d_i_min[4][4][2], double u_d_i_max[4][4][2], double*** u_dx_c, double*** u_dy_c, double u_d_i[4][4][2], double alpha_i[4])
{
  for(unsigned int sol_i = 0; sol_i < 4; sol_i++)
  {
    alpha_i[sol_i] = 1;
    double u_dx = u_dx_c[sol_i][0][1];
    double u_dy = u_dy_c[sol_i][0][2];

    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
      // Sanity checks.
      if(std::abs(u_dx) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[sol_i][vertex_i][0] - u_dx) / u_dx) > 10)
        continue;
      if(std::abs((u_d_i_max[sol_i][vertex_i][0] - u_dx) / u_dx) > 10)
        continue;

      if(std::abs(u_dy) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[sol_i][vertex_i][1] - u_dy) / u_dy) > 10)
        continue;
      if(std::abs((u_d_i_max[sol_i][vertex_i][1] - u_dy) / u_dy) > 10)
        continue;

      // dx.
      if(u_d_i[sol_i][vertex_i][0] < u_dx)
      {
        if((u_d_i_min[sol_i][vertex_i][0] - u_dx) / (u_d_i[sol_i][vertex_i][0] - u_dx) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_min[sol_i][vertex_i][0] - u_dx) / (u_d_i[sol_i][vertex_i][0] - u_dx);
      }
      else
      {
        if((u_d_i_max[sol_i][vertex_i][0] - u_dx) / (u_d_i[sol_i][vertex_i][0] - u_dx) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_max[sol_i][vertex_i][0] - u_dx) / (u_d_i[sol_i][vertex_i][0] - u_dx);
      }
      // dy.
      if(u_d_i[sol_i][vertex_i][1] < u_dy)
      {
        if((u_d_i_min[sol_i][vertex_i][1] - u_dy) / (u_d_i[sol_i][vertex_i][1] - u_dy) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_min[sol_i][vertex_i][1] - u_dy) / (u_d_i[sol_i][vertex_i][1] - u_dy);
      }
      else
      {
        if((u_d_i_max[sol_i][vertex_i][1] - u_dy) / (u_d_i[sol_i][vertex_i][1] - u_dy) < alpha_i[sol_i])
          alpha_i[sol_i] = (u_d_i_max[sol_i][vertex_i][1] - u_dy) / (u_d_i[sol_i][vertex_i] [1]- u_dy);
      }
    }
  }
}

void KuzminDiscontinuityDetector::find_alpha_i_second_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4][2], double u_dx_c[4], double u_dy_c[4], double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4], double alpha_i_real[4])
{
  double c_x, c_y;
  if(e->get_nvert() == 3)
  {
    c_x = (1/3) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x);
    c_y = (1/3) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y);
  }
  else
  {
    c_x = (1/4) * (e->vn[0]->x + e->vn[1]->x + e->vn[2]->x + e->vn[3]->x);
    c_y = (1/4) * (e->vn[0]->y + e->vn[1]->y + e->vn[2]->y + e->vn[3]->y);
  }

  alpha_i_real[0] = alpha_i_real[1] = alpha_i_real[2] = alpha_i_real[3] = -std::numeric_limits<double>::infinity();
  for(unsigned int sol_i = 0; sol_i < this->solutions.size(); sol_i++)
    for(unsigned int vertex_i = 0; vertex_i < e->get_nvert(); vertex_i++)
    {
      // dxx + dxy.
      if( (u_i[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_dxx_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dxy_c[sol_i] * (e->vn[vertex_i]->y - c_y)) > alpha_i_real[sol_i])
        alpha_i_real[sol_i] = (u_i[sol_i][vertex_i][0] - u_dx_c[sol_i]) / (u_dxx_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dxy_c[sol_i] * (e->vn[vertex_i]->y - c_y));
      // dyy + dxy.
      if( (u_i[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_dyy_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dxy_c[sol_i] * (e->vn[vertex_i]->y - c_y)) > alpha_i_real[sol_i])
        alpha_i_real[sol_i] = (u_i[sol_i][vertex_i][1] - u_dy_c[sol_i]) / (u_dyy_c[sol_i] * (e->vn[vertex_i]->x - c_x) + u_dxy_c[sol_i] * (e->vn[vertex_i]->y - c_y));
    }
}


FluxLimiter::FluxLimiter(FluxLimiter::LimitingType type, double* solution_vector, Hermes::vector<SpaceSharedPtr<double>  > spaces, bool Kuzmin_limit_all_orders_independently) : solution_vector(solution_vector), spaces(spaces), limitOscillations(false)
{
  for(unsigned int sol_i = 0; sol_i < spaces.size(); sol_i++)
    limited_solutions.push_back(new Hermes::Hermes2D::Solution<double>(spaces[sol_i]->get_mesh()));

  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);
  switch(type)
  {
  case Krivodonova:
    this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
    break;
  case Kuzmin:
    this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
    break;
  }
};

FluxLimiter::FluxLimiter(FluxLimiter::LimitingType type, Hermes::vector<MeshFunctionSharedPtr<double> > solutions, Hermes::vector<SpaceSharedPtr<double>  > spaces, bool Kuzmin_limit_all_orders_independently) : spaces(spaces), limitOscillations(false)
{
  for(unsigned int sol_i = 0; sol_i < spaces.size(); sol_i++)
    limited_solutions.push_back(new Hermes::Hermes2D::Solution<double>(spaces[sol_i]->get_mesh()));

  this->solution_vector = new double[Space<double>::get_num_dofs(spaces)];
  OGProjection<double> ogProj;

  ogProj.project_global(this->spaces, solutions, this->solution_vector);
  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);

  switch(type)
  {
  case Krivodonova:
    this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
    break;
  case Kuzmin:
    this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
    break;
  }
};

FluxLimiter::~FluxLimiter()
{
	delete detector;
};

void FluxLimiter::get_limited_solutions(Hermes::vector<MeshFunctionSharedPtr<double> > solutions_to_limit)
{
  for(unsigned int i = 0; i < solutions_to_limit.size(); i++)
    solutions_to_limit[i]->copy(this->limited_solutions[i]);
}

int FluxLimiter::limit_according_to_detector(Hermes::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit)
{
  std::set<int> discontinuous_elements = this->detector->get_discontinuous_element_ids();
  std::set<std::pair<int, double> > oscillatory_element_idsRho = this->detector->get_oscillatory_element_idsRho();
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVX = this->detector->get_oscillatory_element_idsRhoVX();
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVY = this->detector->get_oscillatory_element_idsRhoVY();
  std::set<std::pair<int, double> > oscillatory_element_idsRhoE = this->detector->get_oscillatory_element_idsRhoE();

  // First adjust the solution_vector.
  int running_dofs = 0;
  for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
  {
    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(*it);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
    }
    if(this->limitOscillations)
      running_dofs += spaces[space_i]->get_num_dofs();
  }
  if(limitOscillations)
  {
    unsigned int space_i = 0;
    running_dofs = 0;
    for(std::set<std::pair<int, double> >::iterator it = oscillatory_element_idsRho.begin(); it != oscillatory_element_idsRho.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(it->first);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(it->first), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
        else if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = it->second;
    }
    running_dofs += spaces[space_i]->get_num_dofs();
    space_i = 1;
    for(std::set<std::pair<int, double> >::iterator it = oscillatory_element_idsRhoVX.begin(); it != oscillatory_element_idsRhoVX.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(it->first);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(it->first), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
        else if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = it->second;
    }
    running_dofs += spaces[space_i]->get_num_dofs();
    space_i = 2;
    for(std::set<std::pair<int, double> >::iterator it = oscillatory_element_idsRhoVY.begin(); it != oscillatory_element_idsRhoVY.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(it->first);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(it->first), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
        else if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = it->second;
    }
    running_dofs += spaces[space_i]->get_num_dofs();
    space_i = 3;
    for(std::set<std::pair<int, double> >::iterator it = oscillatory_element_idsRhoE.begin(); it != oscillatory_element_idsRhoE.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(it->first);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(it->first), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
        else if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) == 0)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = it->second;
    }
  }

  // Now adjust the solutions.
  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);

  if(coarse_spaces_to_limit != Hermes::vector<SpaceSharedPtr<double> >()) 
  {
    // Now set the element order to zero.
    Element* e;

    for_all_elements(e, spaces[0]->get_mesh())
      e->visited = false;

    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) 
    {
      e = spaces[0]->get_mesh()->get_element(*it);
      AsmList<double> al;
      spaces[0]->get_element_assembly_list(e, &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++) {
        if(H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0 || H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 0) {
          spaces[0]->get_mesh()->get_element(*it)->visited = true;
          bool all_sons_visited = true;
          for(unsigned int son_i = 0; son_i < 4; son_i++)
            if(!spaces[0]->get_mesh()->get_element(*it)->parent->sons[son_i]->visited)
            {
              all_sons_visited = false;
              break;
            }
            if(all_sons_visited)
              for(unsigned int space_i = 0; space_i < spaces.size(); space_i++) 
                coarse_spaces_to_limit[space_i]->set_element_order(spaces[space_i]->get_mesh()->get_element(*it)->parent->id, 0);
        }
      }
    }

    for(int i = 0; i < coarse_spaces_to_limit.size(); i++)
      coarse_spaces_to_limit.at(i)->assign_dofs();
  }

  return discontinuous_elements.size() + oscillatory_element_idsRho.size() + oscillatory_element_idsRhoVX.size() + oscillatory_element_idsRhoVY.size() + oscillatory_element_idsRhoE.size();
};

void FluxLimiter::limit_second_orders_according_to_detector(Hermes::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit)
{
  std::set<int> discontinuous_elements;
  if(dynamic_cast<KuzminDiscontinuityDetector*>(this->detector))
    discontinuous_elements = static_cast<KuzminDiscontinuityDetector*>(this->detector)->get_second_order_discontinuous_element_ids();
  else
    throw Hermes::Exceptions::Exception("limit_second_orders_according_to_detector() is to be used only with Kuzmin's vertex based detector.");

  // First adjust the solution_vector.
  int running_dofs = 0;
  for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
  {
    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) 
    {
      Element* e = spaces[space_i]->get_mesh()->get_element(*it);
      AsmList<double> al;
      spaces[space_i]->get_element_assembly_list(spaces[space_i]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
        if(H2D_GET_H_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 1 || H2D_GET_V_ORDER(spaces[space_i]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 1)
          solution_vector[running_dofs + al.get_dof()[shape_i]] = 0.0;
    }
    //running_dofs += spaces[space_i]->get_num_dofs();
  }
  // Now adjust the solutions.
  Solution<double>::vector_to_solutions(solution_vector, spaces, limited_solutions);
  if(dynamic_cast<KuzminDiscontinuityDetector*>(this->detector))
  {
    bool Kuzmin_limit_all_orders_independently = dynamic_cast<KuzminDiscontinuityDetector*>(this->detector)->get_limit_all_orders_independently();
    delete detector;
    this->detector = new KuzminDiscontinuityDetector(spaces, limited_solutions, Kuzmin_limit_all_orders_independently);
  }
  else
  {
    delete detector;
    this->detector = new KrivodonovaDiscontinuityDetector(spaces, limited_solutions);
  }

  if(coarse_spaces_to_limit != Hermes::vector<SpaceSharedPtr<double> >()) {
    // Now set the element order to zero.
    Element* e;

    for_all_elements(e, spaces[0]->get_mesh())
      e->visited = false;

    for(std::set<int>::iterator it = discontinuous_elements.begin(); it != discontinuous_elements.end(); it++) {
      AsmList<double> al;
      spaces[0]->get_element_assembly_list(spaces[0]->get_mesh()->get_element(*it), &al);
      for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++) {
        if(H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 1 || H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())) > 1) {
          int h_order_to_set = std::min(1, H2D_GET_H_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())));
          int v_order_to_set = std::min(1, H2D_GET_V_ORDER(spaces[0]->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode())));
          spaces[0]->get_mesh()->get_element(*it)->visited = true;
          bool all_sons_visited = true;
          for(unsigned int son_i = 0; son_i < 4; son_i++)
            if(!spaces[0]->get_mesh()->get_element(*it)->parent->sons[son_i]->visited)
            {
              all_sons_visited = false;
              break;
            }
            if(all_sons_visited)
              for(unsigned int space_i = 0; space_i < spaces.size(); space_i++) 
                coarse_spaces_to_limit[space_i]->set_element_order(spaces[space_i]->get_mesh()->get_element(*it)->parent->id, H2D_MAKE_QUAD_ORDER(h_order_to_set, v_order_to_set));
        }
      }
    }

    for(int i = 0; i < coarse_spaces_to_limit.size(); i++)
      coarse_spaces_to_limit.at(i)->assign_dofs();
  }
};

void MachNumberFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(kappa * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / values.at(0)[i]);
}

void PressureFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (kappa - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}

void VelocityFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = values.at(1)[i] / values.at(0)[i];
}

void EntropyFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / p_ext)
      / Hermes::pow((values.at(0)[i] / rho_ext), kappa));
}
