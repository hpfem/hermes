#include "definitions.h"

//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity: du/dt = div[lambda(u) grad u] + f.
//
//  Nonlinearity: lambda(u) = 1 + pow(u, alpha).
//
//  BC: (Non)constant Dirichlet.
//
//  IC: Custom initial condition matching the BC.
//
    
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;

#pragma region Time-dependency setup.
  // Time step. 
  double time_step_length = 0.05;                           
  // Time interval length.
  double current_time = 0.;
  int time_step_number = 1;
  const double T_FINAL = 10.0;
#pragma endregion

#pragma region Adaptivity setup.
  const CandList CAND_LIST = H2D_HP_ANISO;          
  const double ERR_STOP = .5;
  double error_estimate;
#pragma endregion
  
#pragma region Custom problem attributes initialization.
  // Mesh.
  MeshSharedPtr mesh(new Mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_zero("Zero", 0.);
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> essential_bcs(Hermes::vector<EssentialBoundaryCondition<double>*>(&bc_zero, &bc_essential));

  // Initialize the initial condition.
  MeshFunctionSharedPtr<double> sln_time_prev(new CustomInitialCondition(mesh));

  // Thermal conductivity.
#ifdef LINEAR_NONLINEAR_SWITCH
  Hermes1DFunction<double> lambda(10000.);
#else
  CustomNonlinearity lambda(4.0);
#endif

  // Heat source.
  Hermes2DFunction<double> heat_src(1.0);
    
  // Weak formulation initialization.
  // - steady state
  CustomWeakFormSteadyState weak_formulation(&lambda, &heat_src);
  // - time dependent
  // CustomWeakFormTimeDependent weak_formulation(&lambda, &heat_src, sln_time_prev);
#pragma endregion

#pragma region Hermes attributes initialization.
  // Solutions (fine mesh, coarse mesh, previous time level).
  MeshFunctionSharedPtr<double> sln_time_new(new Solution<double>()), sln_time_new_coarse(new Solution<double>());

  // Reference problem creators.
  Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
  Space<double>::ReferenceSpaceCreator ref_space_creator;

  // Linear - Newton solver.
#ifdef LINEAR_NONLINEAR_SWITCH
  LinearSolver<double> solver;
#else
  NewtonSolver<double> solver;
#endif
  double* Newton_initial_guess = NULL;

  // Error calculation.
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
      
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<double> stoppingCriterion(0.75);

  // Adaptivity processor class.
  Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);

  // Element refinement type refinement_selector.
  H1ProjBasedSelector<double> refinement_selector(CAND_LIST);

  // OpenGL views.
  OrderView order_view("Orders", new WinGeom(835, 0, 520, 440));
  ScalarView solution_view("Solution", new WinGeom(560, 270, 520, 440));
#pragma endregion

int main(int argc, char* argv[])
{
  // Load mesh.
  load_mesh(mesh, "domain.xml", INIT_REF_NUM);
  
  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &essential_bcs, P_INIT));
  solver.set_weak_formulation(&weak_formulation);
  solver.set_space(space);

  #pragma region Time stepping loop.
  /*
    solver.set_time_step(time_step_length);
    do
    {
      std::cout << "Time step: " << time_step_number << std::endl;

      #pragma region Spatial adaptivity loop.
        adaptivity.set_space(space);
        int adaptivity_step = 1;
        do
        {
          std::cout << "Adaptivity step: " << adaptivity_step << std::endl;
          #pragma region Construct globally refined reference mesh and setup reference space.
            MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
            SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space(space, ref_mesh);
            solver.set_space(ref_space);
          #pragma endregion
          try
          {
            // Solving.
            solver.solve(get_initial_Newton_guess(adaptivity_step, &weak_formulation, space, ref_space, sln_time_prev));
            Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, sln_time_new);
          }
          catch(Exceptions::Exception& e)
          {
            std::cout << e.what();
          }

          // Project the fine mesh solution onto the coarse mesh.
          OGProjection<double>::project_global(space, sln_time_new, sln_time_new_coarse); 

          // Calculate element errors and error estimate.
          errorCalculator.calculate_errors(sln_time_new_coarse, sln_time_new);
          double error_estimate = errorCalculator.get_total_error_squared() * 100;
          std::cout << "Error estimate: " << error_estimate << "%" << std::endl;

          // Visualize the solution and mesh.
          display(sln_time_new, ref_space);

          // If err_est too large, adapt the mesh.
          if (error_estimate < ERR_STOP)
            break;
          else 
            adaptivity.adapt(&refinement_selector);
      
          adaptivity_step++;
        }
        while(true);
      #pragma endregion

      #pragma region No adaptivity in space.
        try
        {
          // Solving.
          solver.solve(sln_time_prev);

          // Get the solution for visualization etc. from the coefficient vector.
          Solution<double>::vector_to_solution(solver.get_sln_vector(), space, sln_time_new);
      
          // Visualize the solution and mesh.
          display(sln_time_new, space);
        }
        catch(Exceptions::Exception& e)
        {
          std::cout << e.what();
        }
      #pragma endregion

      sln_time_prev->copy(sln_time_new);

      // Increase current time and counter of time steps.
      current_time += time_step_length;
      time_step_number++;
    }
    while (current_time < T_FINAL);
  */
  #pragma endregion

  #pragma region No time stepping (= stationary problem).
  try
  {
    // Solving.
    solver.solve();

    // Get the solution for visualization etc. from the coefficient vector.
    Solution<double>::vector_to_solution(solver.get_sln_vector(), space, sln_time_new);
      
    // Visualize the solution and mesh.
    display(sln_time_new, space);
  }
  catch(Exceptions::Exception& e)
  {
    std::cout << e.what();
  }
  #pragma endregion

  View::wait();
  return 0;
}


#pragma region Helper functions

  void load_mesh(MeshSharedPtr& mesh, const char* filename, int num_initial_refinements)
  {
    MeshReaderH2DXML mloader;
    mloader.load(filename, mesh);

    // Perform initial mesh refinements.
    for(int i = 0; i < INIT_REF_NUM; i++)
      mesh->refine_all_elements(0, true);
  }

  void display(MeshFunctionSharedPtr<double>& sln_to_display, SpaceSharedPtr<double>& space_to_display)
  {
    char title[100];
    sprintf(title, "Mesh, time %g", current_time);
    order_view.set_title(title);
    order_view.show(space_to_display);

    sprintf(title, "Solution, time %g", current_time);
    solution_view.set_title(title);
    solution_view.show_mesh(false);
    solution_view.show(sln_to_display);
  }

  double* get_initial_Newton_guess(int adaptivity_step, WeakForm<double>* weak_formulation, SpaceSharedPtr<double> space,
                                   SpaceSharedPtr<double> ref_space, MeshFunctionSharedPtr<double> sln)
{
  if(Newton_initial_guess)
    delete [] Newton_initial_guess;

  Newton_initial_guess = new double [ref_space->get_num_dofs()];

  if (time_step_number == 1 && adaptivity_step == 1)
  {
    // In the first step, project the coarse mesh solution.
    NewtonSolver<double> newton_coarse;
    newton_coarse.set_space(space);
    newton_coarse.set_weak_formulation(weak_formulation);
  }

  OGProjection<double>::project_global(ref_space, sln, Newton_initial_guess);
  return Newton_initial_guess;
}
#pragma endregion