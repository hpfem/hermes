#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Initial number of uniform refinements.
const int INIT_REF_NUM = 0;
const int INIT_REF_NUM_BDY = 2;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.2;
// This is a stopping criterion for Adaptivity.
class CustomAdaptStoppingCriterion : public AdaptivityStoppingCriterion<double>
{
public:
  /// Decide if the refinement at hand will be carried out.
  bool add_refinement(ErrorCalculator<double>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i)
  {
    ErrorCalculator<double>::ElementReference& element_reference = error_calculator->element_references[element_inspected_i];
    if(processed_error_squared < (THRESHOLD*THRESHOLD) * error_calculator->errors_squared_sum)
      return true;
    else
    {
      ErrorCalculator<double>::ElementReference& element_reference_previous = error_calculator->element_references[element_inspected_i - 1];
      if(*(element_reference.error) / *(element_reference_previous.error) > 0.9)
        return true;
    }
    return false;
  }
} stoppingCriterion;

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Error measurement type.
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

const double AMPLITUDE = 1.00;
const double FREQUENCY = 800;
const double TIME_STEPS_PER_FREQUENCY = 50;

double time_step = 1 / (FREQUENCY * TIME_STEPS_PER_FREQUENCY);

// Time step.
void calculate_time_step(int iteration)
{
  return;
}

// Derefinement time_step
const int derefinement_period = 250 * TIME_STEPS_PER_FREQUENCY / FREQUENCY;

// Stopping criterion for adaptivity.
double REFINEMENT_THRESHOLD(int iteration)
{
  double max = 30.0;
  double min = 1.0;

  if((iteration > TIME_STEPS_PER_FREQUENCY - 2) && (iteration < TIME_STEPS_PER_FREQUENCY))
  {
    std::cout << "\tCurrent refinement threshold: " << max << std::endl;
    return max;
  }

  double ratio = std::min(1., (iteration / (TIME_STEPS_PER_FREQUENCY)));
  double toReturn = max + (min - max) * ratio;
  
  std::cout << "\tCurrent refinement threshold: " << toReturn << std::endl;
  return toReturn;
}

Hermes::Hermes2D::Views::ScalarView viewSR("Ref. Solution", new WinGeom(10, 320, 500, 300));
Views::OrderView oviewR("Ref. Polynomial orders", new WinGeom(520, 320, 500, 300));

Hermes::Hermes2D::Views::ScalarView viewErrorValue("Error in Value", new WinGeom(10, 10, 500, 300));
Hermes::Hermes2D::Views::ScalarView viewErrorDerivative("Error in Derivative", new WinGeom(520, 10, 500, 300));

void output(Hermes::vector<MeshFunctionSharedPtr<double> > slns, SpaceSharedPtr<double> space_value, Hermes::vector<MeshFunctionSharedPtr<double> > rslns, SpaceSharedPtr<double> rspace_value)
{
  //viewS.set_min_max_range(-AMPLITUDE, AMPLITUDE);
  // viewS.show(slns[0]);
  viewSR.show(rslns[1]);
  //oview.show(space_value);
  oviewR.show(rspace_value);
  // oview.wait_for_keypress();

  // Add entry to DOF and CPU convergence graphs.
  /*
  Linearizer lin;
  Orderizer ord;

  char* filename = new char[1000];
  sprintf(filename, "Solution-%i.vtk", iteration);
  lin.save_solution_vtk(rsln_value, filename, "sln", true);
  sprintf(filename, "Orders-%i.vtk", iteration);
  ord.save_orders_vtk(rspace_value, filename);
  sprintf(filename, "Mesh-%i.vtk", iteration);
  ord.save_mesh_vtk(rspace_value, filename);
  */
}

class MySelector : public H1ProjBasedSelector<double>
{
public:
  MySelector(CandList cand_list) : H1ProjBasedSelector<double>(cand_list)
  {
    this->set_error_weights(4.0, 1.0, 16.0);
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

  Hermes::vector<Cand> create_candidates(Element* e, int quad_order)
  {
    Hermes::vector<Cand> candidates;

    if(this->cand_list == H2D_NONE)
    {

      // Get the current order range.
      int current_min_order, current_max_order;
      this->get_current_order_range(e, current_min_order, current_max_order);

      int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

      if(current_max_order < std::max(order_h, order_v))
        current_max_order = std::max(order_h, order_v);

      int last_order_h = std::min(current_max_order, order_h + 1), last_order_v = std::min(current_max_order, order_v + 1);
      int last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);

      candidates.push_back(Cand(H2D_REFINEMENT_P, quad_order));
      candidates.push_back(Cand(H2D_REFINEMENT_P, last_order));
      candidates.push_back(Cand(H2D_REFINEMENT_H, quad_order, quad_order, quad_order, quad_order));
      
      return candidates;
    }
    else
    {
      return H1ProjBasedSelector<double>::create_candidates(e, quad_order);
    }
  }
};

class MyErrorCalculator : public ErrorCalculator<double>
{
  public:
  class MyNormFormVol : public NormFormVol<double>
  {
  public:
    MyNormFormVol(int i, int j) : NormFormVol<double>(i, j) {};
      
    double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
    {
      double result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * ((u->val[i] * conj(v->val[i])) + (u->dx[i] * conj(v->dx[i])) + (u->dy[i] * conj(v->dy[i])));
      return result * e->area;
    }
  };

  MyErrorCalculator(CalculatedErrorType errorType) : ErrorCalculator(errorType)
  {
    for(int i = 0; i < 2; i++)
    {
      this->add_error_form(new MyNormFormVol(i, i));
    }
  }

  virtual ~MyErrorCalculator() {}
};

void load_mesh(MeshSharedPtr mesh_value, MeshSharedPtr mesh_derivative, MeshSharedPtr basemesh)
{
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh_value);
  try
  {
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("test.xml", mesh_value);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  for(int i = 0; i < INIT_REF_NUM; i++)
    mesh_value->refine_all_elements(0, true);

  mesh_value->refine_towards_boundary("Source", INIT_REF_NUM_BDY, false, false);
  mesh_value->refine_towards_vertex(0, INIT_REF_NUM_BDY, false);
  mesh_value->refine_towards_vertex(1, INIT_REF_NUM_BDY, false);
  
  basemesh->copy(mesh_value);
  mesh_derivative->copy(mesh_value);
}
int main(int argc, char* argv[])
{
  //Hermes2DApi.set_integral_param_value(numThreads, 1);

  Hermes::Mixins::Loggable logger;

  logger.set_verbose_output(true);

  // Load the mesh.
  MeshSharedPtr mesh_value(new Mesh), mesh_derivative(new Mesh), basemesh(new Mesh);
  load_mesh(mesh_value, mesh_derivative, basemesh);
  
  std::string matched_boundary("Wall");
  double matched_boundaries_values(345*1.25);
  std::string source_boundary("Source");

  CustomBCValue custom_bc_value(source_boundary, AMPLITUDE, FREQUENCY);
  Hermes::Hermes2D::EssentialBCs<double> bcs_value(&custom_bc_value);

  CustomBCDerivative custom_bc_derivative(source_boundary, AMPLITUDE, FREQUENCY);
  Hermes::Hermes2D::EssentialBCs<double> bcs_derivative(&custom_bc_derivative);

  SpaceSharedPtr<double> space_value(new Hermes::Hermes2D::H1Space<double>(mesh_value, &bcs_value, P_INIT));
  SpaceSharedPtr<double> space_derivative(new Hermes::Hermes2D::H1Space<double>(mesh_derivative, &bcs_derivative, P_INIT));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_value, space_derivative);

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln_value(new ZeroSolution<double>(mesh_value));
  MeshFunctionSharedPtr<double> sln_derivative(new ZeroSolution<double>(mesh_derivative));
  MeshFunctionSharedPtr<double> rsln_value(new Solution<double>()), rsln_derivative(new Solution<double>());
  MeshFunctionSharedPtr<double> prev_value(new ZeroSolution<double>(mesh_value)), prev_derivative(new ZeroSolution<double>(mesh_derivative));
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_value, sln_derivative);
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_value, rsln_derivative);
  Hermes::vector<MeshFunctionSharedPtr<double> > prevslns(prev_value, prev_derivative);

  MyWeakForm wf(matched_boundary, matched_boundaries_values, prevslns);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, spaces);
  //linear_solver.set_do_not_use_cache();
  linear_solver.set_jacobian_constant();
  // linear_solver.set_UMFPACK_output(true, true);

  // Adaptivity.
  MyErrorCalculator error_calculator(errorType);
  Adapt<double> adaptivity(spaces, &error_calculator);
  adaptivity.set_strategy(&stoppingCriterion);
  adaptivity.set_regularization_level(MESH_REGULARITY);

  // Initialize refinement selector.
  MySelector selector(CAND_LIST);
  Hermes::vector<Selector<double>*> selectors(&selector, &selector);

  // Solve the linear problem.
  double time = 0.;
  bool force_derefinement = false;
  for(int iteration = 0; ; iteration++)
  {
    calculate_time_step(iteration);
    time += time_step;
    wf.set_current_time_step(time_step);

    logger.info("\n\nIteration: %u, Time: %g.", iteration, time);

    // Periodic global derefinement.
    if ((iteration > 1 && iteration % derefinement_period == 0) || force_derefinement)
    {
      logger.info("\t! Global mesh derefinement.");
      mesh_value->copy(basemesh);
      mesh_derivative->copy(basemesh);
      space_value->set_uniform_order(P_INIT);
      space_derivative->set_uniform_order(P_INIT);
      Space<double>::assign_dofs(Hermes::vector<SpaceSharedPtr<double> >(space_value, space_derivative));
      Element* e;
      for_all_active_elements(e, mesh_value)
      {
        space_value->edata[e->id].changed_in_last_adaptation = false;
      }
      for_all_active_elements(e, mesh_derivative)
      {
        space_derivative->edata[e->id].changed_in_last_adaptation = false;
      }
      linear_solver.set_report_cache_hits_and_misses(true);
      linear_solver.free_cache();
      force_derefinement = false;
    }

    Space<double>::update_essential_bc_values(spaces, time);
    
    int as = 1;
    while (true)
    {
      logger.info("\n\tAdaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space->
      Mesh::ReferenceMeshCreator refMeshCreatorValue(mesh_value);
      MeshSharedPtr ref_mesh_value = refMeshCreatorValue.create_ref_mesh();
      Mesh::ReferenceMeshCreator refMeshCreatorDerivative(mesh_derivative);
      MeshSharedPtr ref_mesh_derivative = refMeshCreatorDerivative.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorValue(space_value, ref_mesh_value);
      SpaceSharedPtr<double> rspace_value = refSpaceCreatorValue.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorDerivative(space_derivative, ref_mesh_derivative);
      SpaceSharedPtr<double> rspace_derivative = refSpaceCreatorDerivative.create_ref_space();
      Hermes::vector<SpaceSharedPtr<double> > rspaces(rspace_value, rspace_derivative);
      Space<double>::update_essential_bc_values(rspaces, time);
      
      linear_solver.set_spaces(rspaces);
      if((iteration % derefinement_period) && as == 1)
        linear_solver.jacobian_reusable = true;

      try
      {
        linear_solver.solve();
        
        double facSize = linear_solver.get_UMFPACK_reporting_data(Solver<double>::FactorizationSize);
        double memSize = linear_solver.get_UMFPACK_reporting_data(Solver<double>::PeakMemoryUsage);
        double flops = linear_solver.get_UMFPACK_reporting_data(Solver<double>::Flops);
        
        linear_solver.set_report_cache_hits_and_misses(false);
      
        Solution<double>::vector_to_solutions(linear_solver.get_sln_vector(), rspaces, rslns);
        OGProjection<double>::project_global(spaces, rslns, slns);

        output(slns, spaces[0], rslns, rspaces[0]);

        error_calculator.calculate_errors(slns, rslns);

        MeshFunctionSharedPtr<double> error_value(new ExactSolutionConstantArray<double>(mesh_value, error_calculator.errors[0])); 
        MeshFunctionSharedPtr<double> error_derivative(new ExactSolutionConstantArray<double>(mesh_derivative, error_calculator.errors[1])); 
        viewErrorValue.show(error_value);
        viewErrorDerivative.show(error_derivative);
      }
      catch(std::exception& e)
      {
        iteration = iteration - 1;
        time = time - time_step;
        break;
      }

      double err_est_rel = error_calculator.get_total_error_squared()* 100.0;
      // Report results.
      logger.info("\tCoarse DOF: %d,  Fine DOF: %d.", Space<double>::get_num_dofs(spaces), Space<double>::get_num_dofs(rspaces));
      logger.info("\tTotal error: %g%%.", err_est_rel);

      // We are above the error.
      if(err_est_rel < REFINEMENT_THRESHOLD(iteration))
        break;
      
      {
        adaptivity.adapt(selectors);
        as++;
      }
    }

    prev_value->copy(rsln_value);
    prev_derivative->copy(rsln_derivative);
  }

  return 0;
}