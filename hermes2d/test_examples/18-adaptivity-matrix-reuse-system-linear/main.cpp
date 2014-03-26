#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

const int P_INIT = 1;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;               // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20;        // Fixed temperature on the boundary.

/// Custom selector, selects the appropriate elements
class CustomSelector : public Selector<double>
{
public:
  /// Constructor.
  CustomSelector() : Selector<double>()
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
      refinement.refinement_polynomial_order[0] = refinement.refinement_polynomial_order[1] = refinement.refinement_polynomial_order[2] = refinement.refinement_polynomial_order[3] = P_INIT;
      return true;
    }
    return false;
  }

  std::vector<int> element_ids;
};

#define ADAPT_SOLVER

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(numThreads, 1);
  //HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  //HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_MUMPS);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), mesh1(new Mesh), mesh2(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("quad.xml", mesh);

  // Refine all elements, do it INIT_REF_NUM-times.
  for (unsigned int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  mesh1->copy(mesh);
  mesh2->copy(mesh);

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential0(std::vector<std::string>({ "0", "1"}), FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs0(&bc_essential0);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential1(std::vector<std::string>({ "0", "1"}), FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs1(&bc_essential1);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential2(std::vector<std::string>({ "0", "1", "2" }), FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs2(&bc_essential1);

  // Initialize space->
  SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs0, P_INIT));
  SpaceSharedPtr<double> space1(new Hermes::Hermes2D::H1Space<double>(mesh1, &bcs1, P_INIT));
  SpaceSharedPtr<double> space2(new Hermes::Hermes2D::H1Space<double>(mesh2, &bcs2, P_INIT));

  std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

  // Weak Form
  WeakFormSharedPtr<double> wf(new CustomWeakFormPoisson(new Hermes::Hermes2DFunction<double>(VOLUME_HEAT_SRC)));

#ifdef ADAPT_SOLVER
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(CalculatedErrorType::RelativeErrorToGlobalNorm, 2);
  AdaptStoppingCriterionSingleElement<double> criterion(0.);
  CustomSelector selector[3];
  for (int j = 0; j < 3; j++)
  {
    for (int i = 0; i < 1000; i++)
    {
      if (((rand() % 2) + i) % 5 == 0)
        selector[j].element_ids.push_back(i);
    }
    selector[j].element_ids.push_back(j);
  }
  AdaptSolverCriterionFixed global_criterion(22);

  AdaptSolver<double, LinearSolver<double> > adaptSolver({ space, space1 }, wf, &errorCalculator, &criterion, { &selector[0], &selector[1] }, &global_criterion);

  adaptSolver.switch_visualization(true);
  adaptSolver.set_verbose_output(true);
  adaptSolver.solve(hpAdaptivity);

#else

  LinearSolver<double> solver(wf, { space, space1 });
  //solver.set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
  //solver.output_matrix();
  solver.solve();
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  Solution<double>::vector_to_solution(solver.get_sln_vector(), space, sln);

  ScalarView s;
  s.get_linearizer()->set_criterion(LinearizerCriterionFixed(5));
  s.show(sln);

#endif
  View::wait();
  return 0;
}