#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Initial polynomial degree for u.
const int P_INIT_U = 4;
// Initial polynomial degree for v.
const int P_INIT_V = 2;
const int INIT_REF_NUM = 5;               // Number of initial uniform mesh refinements.
// Number of initial boundary refinements
const int INIT_REF_BDY = 6;

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100.;

/// Custom selector, selects the appropriate elements
class CustomSelector : public Selector<double>
{
public:
  /// Constructor.
  CustomSelector(int space_i) : Selector<double>(), space_i(space_i)
  {
  };

  CustomSelector(std::vector<int> element_ids) : Selector<double>(), element_ids(element_ids)
  {
  };

  virtual bool select_refinement(Element* element, int quad_order, MeshFunction<double>* rsln, ElementToRefine& refinement)
  {
    for (unsigned short i = 0; i < this->element_ids.size(); i++)
    if (element->id == this->element_ids[i])
    {
      refinement.split = H2D_REFINEMENT_H;
      if (space_i == 0)
        refinement.refinement_polynomial_order[0] = refinement.refinement_polynomial_order[1] = refinement.refinement_polynomial_order[2] = refinement.refinement_polynomial_order[3] = P_INIT_U;
      else
        refinement.refinement_polynomial_order[0] = refinement.refinement_polynomial_order[1] = refinement.refinement_polynomial_order[2] = refinement.refinement_polynomial_order[3] = P_INIT_V;
      return true;
    }
    return false;
  }
  int space_i;
  std::vector<int> element_ids;
};

#define ADAPT_SOLVER

int main(int argc, char* argv[])
{
  //HermesCommonApi.set_integral_param_value(numThreads, 1);
  //HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);

  // Load the mesh.
  MeshSharedPtr u_mesh(new Mesh), v_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", u_mesh);
  
  // Refine all elements, do it INIT_REF_NUM-times.
  for (unsigned int i = 0; i < INIT_REF_NUM; i++)
    u_mesh->refine_all_elements();
  u_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  v_mesh->copy(u_mesh);

// Set exact solutions.
  MeshFunctionSharedPtr<double> exact_u(new ExactSolutionFitzHughNagumo1(u_mesh));
  MeshFunctionSharedPtr<double> exact_v(new ExactSolutionFitzHughNagumo2(v_mesh, K));

  // Define right-hand sides.
  CustomRightHandSide1 g1(K, D_u, SIGMA);
  CustomRightHandSide2 g2(K, D_v);

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u("Bdy", 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v("Bdy", 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(u_mesh, &bcs_u, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(v_mesh, &bcs_v, P_INIT_V));

  // Weak Form
  WeakFormSharedPtr<double> wf(new CustomWeakForm(&g1, &g2));

#ifdef ADAPT_SOLVER
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(CalculatedErrorType::RelativeErrorToGlobalNorm, 2);
  AdaptStoppingCriterionSingleElement<double> criterion(0.);
  CustomSelector selector_u(0);
  for (int i = 0; i < u_mesh->get_max_element_id(); i++)
  {
    if ((i % 33) == 0) selector_u.element_ids.push_back(i);
  }
  CustomSelector selector_v(1);
  for (int i = 0; i < v_mesh->get_max_element_id(); i++)
  {
    if ((i % 33) == 0) selector_v.element_ids.push_back(i);
  }
  AdaptSolverCriterionFixed global_criterion(1);

  AdaptSolver<double, NewtonSolver<double> > adaptSolver({ u_space, v_space }, wf, &errorCalculator, &criterion, { &selector_u, &selector_v }, &global_criterion);

  adaptSolver.switch_visualization(false);
  adaptSolver.set_verbose_output(true);

  try
  {
    adaptSolver.get_solver()->get_linear_matrix_solver()->as_IterSolver()->set_tolerance(1e-5, LoopSolverToleranceType::RelativeTolerance);
    adaptSolver.get_solver()->get_linear_matrix_solver()->as_IterSolver()->set_max_iters(400);
    adaptSolver.get_solver()->get_linear_matrix_solver()->as_IterSolver()->set_solver_type(IterSolverType::CG);
    adaptSolver.get_solver()->get_linear_matrix_solver()->as_IterSolver()->set_precond(new ParalutionPrecond<double>(PreconditionerType::MultiColoredILU));
  }
  catch(std::exception& e)
  {
  }
  adaptSolver.solve(hpAdaptivity);

#else

  LinearSolver<double> solver(wf, space);
  //solver.set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
  //solver.output_matrix();
  solver.solve();
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  Solution<double>::vector_to_solution(solver.get_sln_vector(), space, sln);

  ScalarView s;
  s.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  s.show(sln);

  MeshSharedPtr meshf(new Mesh);
  meshf->copy(mesh);
  SpaceSharedPtr<double> spacef(new H1Space<double>(meshf, &bcs, 1));
  spacef->set_element_order(1, 5);
  spacef->set_element_order(2, 2);
  spacef->set_element_order(3, 2);
  spacef->set_element_order(4, 2);
  spacef->set_element_order(5, 2);
  spacef->assign_dofs();

  BaseView<double> bf("Fine space");
  bf.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  bf.show(spacef);

  WeakFormSharedPtr<double> wff(new WeakFormsH1::DefaultWeakFormPoissonLinear<double>(HERMES_ANY, new Hermes2DFunction<double>(13.0)));

  LinearSolver<double> solverf(wff, spacef);
  solverf.set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
  solverf.output_matrix();
  solverf.solve();

  MeshFunctionSharedPtr<double> slnf(new Solution<double>);

  Solution<double>::vector_to_solution(solverf.get_sln_vector(), spacef, slnf);

  ScalarView sf;
  sf.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  sf.show(slnf);

#endif
  View::wait();
  return 0;
}