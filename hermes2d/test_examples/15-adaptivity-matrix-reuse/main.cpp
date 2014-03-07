#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/// Custom selector, selects the appropriate elements
class CustomSelector : public Selector<double>
{
public:
  /// Constructor.
  CustomSelector(std::vector<int> element_ids) : Selector<double>(), element_ids(element_ids)
  {
  };

  virtual bool select_refinement(Element* element, int quad_order, MeshFunction<double>* rsln, ElementToRefine& refinement)
  {
    for(unsigned short i = 0; i < this->element_ids.size(); i++)
    if (element->id == this->element_ids[i])
    {
      refinement.refinement_polynomial_order[0] = refinement.refinement_polynomial_order[1] = refinement.refinement_polynomial_order[2] = refinement.refinement_polynomial_order[3] = 1;
      return true;
    }
    return false;
  }

  std::vector<int> element_ids;
};

#define ADAPT_SOLVER

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("quad.xml", mesh);
  mesh->refine_all_elements();
  mesh->refine_towards_boundary("0", 4);
  mesh->refine_towards_boundary("1", 3);
  mesh->refine_towards_vertex(154,2);

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("0", 10.);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space.
  SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, 1));

  BaseView<double> b("Coarse space");
  b.get_linearizer()->set_criterion(LinearizerCriterionFixed(5));
  b.show(space);

  // Weak Form
  WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormPoissonLinear<double>(HERMES_ANY, new Hermes2DFunction<double>(13.0)));

#ifdef ADAPT_SOLVER
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(CalculatedErrorType::RelativeErrorToGlobalNorm, 1);
  AdaptStoppingCriterionSingleElement<double> criterion(0.);
  CustomSelector selector({ 0, 2, 10, 34, 132, 351, 400, 500, 600, 700, 800 });
  AdaptSolverCriterionFixed global_criterion(2);

  AdaptSolver<double, LinearSolver<double> > adaptSolver(space, wf, &errorCalculator, &criterion, &selector, &global_criterion);

  adaptSolver.switch_visualization(true);
  adaptSolver.set_verbose_output(true);
  adaptSolver.solve();

  ScalarView s;
  s.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  s.show(adaptSolver.get_ref_sln(0));

#else

  LinearSolver<double> solver(wf, space);
  solver.set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
  solver.output_matrix();
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

  View::wait();
#endif
  return 0;
}