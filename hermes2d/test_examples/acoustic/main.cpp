#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Tiem step.
const double time_step = 3e-5;
// Initial number of uniform refinements.
const int INIT_REF_NUM = 2;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.666;
// This is a stopping criterion for Adaptivity.
const AdaptivityStoppingCriterion stoppingCriterion = AdaptStoppingCriterionSingleElement;
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Error measurement type.
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

// Stopping criterion for adaptivity.
double REFINEMENT_THRESHOLD(double time, double norm)
{
  double toReturn = 5.0;
  /*
    if(time > time_step * 10)
    {
      if(std::log10(norm / 1e10) > 0)
        toReturn = toReturn / (std::log10(norm / 1e10) + 1);
    }

  */
  
  toReturn = toReturn / std::max(std::log10(time / time_step) / 3, 1.0);

  std::cout << "\tCurrent refinement threshold: " << toReturn << std::endl;

  return toReturn;
}

double AMPLITUDE = 1.0;
double FREQUENCY = 1000;

class MySelector : public H1ProjBasedSelector<double>
{
public:
  MySelector(CandList cand_list) : H1ProjBasedSelector<double>(cand_list)
  {
  }
private:
  void evaluate_cands_score(Hermes::vector<Cand>& candidates, Hermes2D::Element* e)
  {
    //calculate score of candidates
    Cand& unrefined = candidates[0];
    const int num_cands = (int)candidates.size();
    unrefined.score = 0;

    for (int i = 1; i < num_cands; i++)
    {
      Cand& cand = candidates[i];
      if(cand.error < unrefined.error)
      {
        double delta_dof = cand.dofs - unrefined.dofs;
        candidates[i].score = (log10(unrefined.error) - log10(cand.error)) / delta_dof;
      }
      else
        candidates[i].score = 0;
    }
  }
};

int main(int argc, char* argv[])
{
  //Hermes2DApi.set_integral_param_value(numThreads, 1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);
  try
  {
    // Hermes::Hermes2D::MeshReaderH2D mloader;
    // mloader.load("domain_apartment.mesh", mesh);

    //Hermes::Hermes2D::MeshReaderH2DXML mloader;
    //mloader.load("domain_apartment_test.xml", mesh);

    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("test.xml", mesh);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  /*
     MeshView m;
     m.show(mesh);
     m.wait_for_close();
  */

  for(int i = 0; i < INIT_REF_NUM; i++)
  {
    mesh->refine_all_elements(0, true);
    mesh->refine_towards_boundary("Source", 1, true, true);
  }
  
  std::string matched_boundary("Wall");
  double matched_boundaries_values(345*1.25);
  std::string source_boundary("Source");

  CustomBCValue custom_bc_value(source_boundary, AMPLITUDE, FREQUENCY);
  Hermes::Hermes2D::EssentialBCs<double> bcs_value(&custom_bc_value);

  CustomBCDerivative custom_bc_derivative(source_boundary, AMPLITUDE, FREQUENCY);
  Hermes::Hermes2D::EssentialBCs<double> bcs_derivative(&custom_bc_derivative);

  SpaceSharedPtr<double> space_value(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs_value, P_INIT));
  SpaceSharedPtr<double> space_derivative(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs_derivative, P_INIT));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_value, space_derivative);

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln_value(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_derivative(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> rsln_value(new Solution<double>()), rsln_derivative(new Solution<double>());
  MeshFunctionSharedPtr<double> prev_value(new ZeroSolution<double>(mesh)), prev_derivative(new ZeroSolution<double>(mesh));
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_value, sln_derivative);
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_value, rsln_derivative);
  Hermes::vector<MeshFunctionSharedPtr<double> > prevslns(prev_value, prev_derivative);

  Hermes::Hermes2D::Views::ScalarView viewS("Solution");
  viewS.set_min_max_range(-AMPLITUDE, AMPLITUDE);
  Views::OrderView oview("Polynomial orders");

  MyWeakForm wf(matched_boundary, matched_boundaries_values, prevslns);
  wf.set_current_time_step(time_step);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, spaces);
  linear_solver.set_do_not_use_cache();
  linear_solver.set_jacobian_constant();

  // Adaptivity.
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 2);
  Adapt<double> adaptivity(spaces, &error_calculator);
  adaptivity.set_strategy(stoppingCriterion, THRESHOLD);
  // Initialize refinement selector.
  MySelector selector(CAND_LIST);
  Hermes::vector<Selector<double>*> selectors(&selector, &selector);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est;
  
  // Solve the linear problem.
  unsigned int iteration = 0;

  // Unrefine when derefinement_threshold is reached.
  bool derefinement_threshold_met = false;

  for(double time = time_step; time < 1.0; time += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("\n\nIteration: %u, Time: %g.", iteration, time);
    Space<double>::update_essential_bc_values(spaces, time);

    // Periodic global derefinement.
    bool derefined_in_this_step = false;
    if (derefinement_threshold_met)
    {
      Hermes::Mixins::Loggable::Static::info("\t! Global mesh derefinement.");
      mesh->unrefine_all_elements();
      space_value->adjust_element_order(-1, -1, P_INIT, P_INIT);
      space_derivative->adjust_element_order(-1, -1, P_INIT, P_INIT);
      Space<double>::assign_dofs(Hermes::vector<SpaceSharedPtr<double> >(space_value, space_derivative));
      derefinement_threshold_met = false;
      derefined_in_this_step = true;
    }

    int as = 1;
    while (true)
    {
      // Construct globally refined reference mesh and setup reference space->
      Mesh::ReferenceMeshCreator refMeshCreator(mesh);
      MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreator(space_value, ref_mesh);
      SpaceSharedPtr<double> rspace_value = refSpaceCreator.create_ref_space();
      SpaceSharedPtr<double> rspace_derivative = refSpaceCreator.create_ref_space();
      Hermes::vector<SpaceSharedPtr<double> > rspaces(rspace_value, rspace_derivative);
      Space<double>::update_essential_bc_values(rspaces, time);
      
      Hermes::Mixins::Loggable::Static::info("\n\tAdaptivity step %d:", as++);
      Hermes::Mixins::Loggable::Static::info("\tCoarse DOF: %d,  Fine DOF: %d.", Space<double>::get_num_dofs(spaces), Space<double>::get_num_dofs(rspaces));
      
      linear_solver.set_spaces(rspaces);
      linear_solver.solve();
      
      Solution<double>::vector_to_solutions(linear_solver.get_sln_vector(), rspaces, rslns);
      OGProjection<double>::project_global(spaces, rslns, slns);
      error_calculator.calculate_errors(slns, rslns);
      double err_est_rel = error_calculator.get_error_squared(0) * 100.0;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("\tTotal error: %g%%.", err_est_rel);
      Hermes::Mixins::Loggable::Static::info("\tTotal norm: %g.", error_calculator.get_total_norm_squared());

      // View the coarse mesh solution and polynomial orders.
      viewS.show(rsln_value);
      oview.show(rspace_value);

      // Add entry to DOF and CPU convergence graphs.
      /*
      graph_dof_est.add_values(rspace_value->get_num_dofs(), err_est_rel);
      graph_dof_est.save("conv_dof_est.dat");

      Linearizer lin;
      char* filename = new char[1000];
      sprintf(filename, "Solution-%i.vtk", iteration);
      lin.save_solution_vtk(rsln_value, filename, "sln", true);
      Orderizer ord;
      sprintf(filename, "Orders-%i.vtk", iteration);
      ord.save_orders_vtk(rspace_value, filename);
      sprintf(filename, "Mesh-%i.vtk", iteration);
      ord.save_mesh_vtk(rspace_value, filename);
      */

      if(derefinement_threshold_met && err_est_rel < REFINEMENT_THRESHOLD(time, error_calculator.get_total_norm_squared()))
        break;

      while(err_est_rel < REFINEMENT_THRESHOLD(time, error_calculator.get_total_norm_squared()))
      {
        prev_value->copy(rsln_value);
        prev_derivative->copy(rsln_derivative);
        time += time_step;
        derefined_in_this_step = false;
        as = 1;
        Hermes::Mixins::Loggable::Static::info("\n\nIteration - inner loop: %u, Time: %g.", ++iteration, time);
        Hermes::Mixins::Loggable::Static::info("\tCoarse DOF: %d,  Fine DOF: %d.", Space<double>::get_num_dofs(spaces), Space<double>::get_num_dofs(rspaces));
        Space<double>::update_essential_bc_values(rspaces, time);
        
        linear_solver.solve();
        
        Solution<double>::vector_to_solutions(linear_solver.get_sln_vector(), rspaces, rslns);
        OGProjection<double>::project_global(spaces, rslns, slns);
        error_calculator.calculate_errors(slns, rslns);
        err_est_rel = error_calculator.get_error_squared(0) * 100.0;
        
        // Report results.
        Hermes::Mixins::Loggable::Static::info("\tTotal error: %g%%.", err_est_rel);
        Hermes::Mixins::Loggable::Static::info("\tTotal norm: %g.", error_calculator.get_total_norm_squared());

        viewS.show(rsln_value);
        oview.show(rspace_value);
      }
      
      // We are above the error.
      {
        adaptivity.adapt(selectors);
        if(!derefined_in_this_step)
          derefinement_threshold_met = true;
      }
    }

    prev_value->copy(rsln_value);
    prev_derivative->copy(rsln_derivative);
    iteration++;
  }

  return 0;
}