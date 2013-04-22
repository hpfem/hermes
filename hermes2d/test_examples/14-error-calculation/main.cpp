#include "definitions.h"

#define CUSTOM_DEBUG

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh_coarse(new Mesh), mesh_fine(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh_coarse);

  mesh_fine->copy(mesh_coarse);
  mesh_fine->refine_all_elements();

  DefaultEssentialBCConst<double> bc_coarse("Bdy", 1.0);
  DefaultEssentialBCConst<double> bc_fine("Bdy", 2.0);

  EssentialBCs<double> bcs_coarse(&bc_coarse);
  EssentialBCs<double> bcs_fine(&bc_fine);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space_coarse(new H1Space<double>(mesh_coarse, &bcs_coarse, 1));
  SpaceSharedPtr<double> space_fine(new H1Space<double>(mesh_fine, &bcs_fine, 1));

  // Translate the coefficient vector into a Solution. 
  MeshFunctionSharedPtr<double> fn_coarse(new ConstantSolution<double>(mesh_coarse, 1.0));
  MeshFunctionSharedPtr<double> fn_fine_const(new ConstantSolution<double>(mesh_fine, 2.0));
  MeshFunctionSharedPtr<double> fn_fine_non_const(new CustomExactSolutionScalar(mesh_fine));
 
#ifdef CUSTOM_DEBUG
  // Visualise the solution and mesh.
  ScalarView s_view;
  s_view.show(fn_coarse);
  s_view.wait_for_keypress();
  s_view.show(fn_fine_const);
  s_view.wait_for_keypress();
  s_view.show(fn_fine_non_const);
  s_view.wait_for_keypress();
#endif

  ErrorCalculator<double> errorCalculator(AbsoluteError);
  errorCalculator.add_error_form(new CustomNormFormVol(0,0));
  errorCalculator.add_error_form(new CustomNormFormSurf(0,0));
  errorCalculator.add_error_form(new CustomNormFormDG(0,0));

  errorCalculator.calculate_errors(fn_coarse, fn_fine_const);

#ifdef CUSTOM_DEBUG
  std::cout << "Total error const: " << errorCalculator.get_total_error_squared() << std::endl;
#endif

  if(std::abs(errorCalculator.get_total_error_squared() - 5.0) > 1e-10)
    return -1;

  errorCalculator.calculate_errors(fn_coarse, fn_fine_non_const);

#ifdef CUSTOM_DEBUG
  std::cout << "Total error nonconst: " << errorCalculator.get_total_error_squared() << std::endl;
#endif

  if(std::abs(errorCalculator.get_total_error_squared() - 13.0) > 1e-10)
    return -1;

  std::cout << "Success!";
  return 0;
}