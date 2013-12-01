#include "definitions.h"
#include "../euler_util.h"
#include "algorithms.h"

int polynomialDegree = 2;
int initialRefinementsCount = 4;
const Algorithm algorithm = Multiscale;
const SolvedExample solvedExample = SolvedExample::MovingPeak;
double MovingPeakDiffusivity = 1e-2;
const EulerLimiterType limiter_type = VertexBased;

double diffusivity = 1e-3;
double s = -1;
double CFL = 1.;

int main(int argc, char* argv[])
{
  if(argc > 1)
    initialRefinementsCount = atoi(argv[1]);
  if(argc > 2)
    diffusivity = (double)atof(argv[2]);
  if(argc > 3)
    CFL = atof(argv[3]);
    
  double sigma = std::pow(2., (double)(initialRefinementsCount)) * (s == -1 ? 10.0 : (s == 1 ? 10. : 0.));

  // test();
  Hermes::Mixins::Loggable logger(true);
  Hermes::Mixins::Loggable logger_details(true);
  std::stringstream ss;
  ss << "logfile_p=" << polynomialDegree << "_" << initialRefinementsCount << "_eps=" << diffusivity << "_s=" << s << ".h2d";
  logger.set_logFile_name(ss.str());
  std::stringstream ssd;
  ssd << "logfile_detail_p=" << polynomialDegree << "_" << initialRefinementsCount << "_eps=" << diffusivity << "_s=" << s << ".h2d";
  logger_details.set_logFile_name(ssd.str());
  
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  double mesh_size;

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    mloader.load("domain_rotation.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
    mesh_size = 1.;
  case AdvectedCube:
    mloader.load("larger_domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 1.;
    break;
  case CircularConvection:
    mloader.load("domain_circular_convection.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 1.;
    break;
  case MovingPeak:
    mloader.load("domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = 2.;
    break;
  case Benchmark:
    mloader.load("domain_benchmark.xml", mesh);
    mesh->refine_all_elements(2);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    mesh_size = .6;
    break;
  }

  double time_step_length = mesh_size * CFL * std::pow(2., -(double)initialRefinementsCount);
  logger.info("Time step: %f.", time_step_length);

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition_der = NULL;
  ExactSolutionScalar<double>* initial_solution = NULL;
  ExactSolutionScalar<double>* exact_sln = NULL;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
    initial_condition_der = new InitialConditionSolidBodyRotation(mesh);
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
    break;
  case AdvectedCube:
    initial_condition = new InitialConditionAdvectedCube(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);
    break;
  case CircularConvection:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    exact_sln = new ExactSolutionCircularConvection(mesh);
    initial_solution = new ExactSolutionCircularConvection(mesh);
    break;
  case MovingPeak:
    initial_condition = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    initial_condition_der = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    previous_initial_condition = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.);
    exact_sln = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, (1./2.) * M_PI);
    initial_solution = new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, (1./2.) * M_PI);
    break;
  case Benchmark:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    exact_sln = new ExactSolutionBenchmark2(mesh, diffusivity);
    initial_solution = new ExactSolutionBenchmark2(mesh, diffusivity);
    break;
  }
  
  // Solutions.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);
  MeshFunctionSharedPtr<double> previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double> previous_mean_values(initial_condition);
  MeshFunctionSharedPtr<double> previous_derivatives(initial_condition_der);
  MeshFunctionSharedPtr<double> exact_solution(exact_sln);
  MeshFunctionSharedPtr<double> initial_sln(initial_solution);

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(0, 0, 600, 350));
  ScalarView exact_view("Exact solution", new WinGeom(610, 0, 600, 350));
  exact_view.show(exact_solution);
  
  // Exact solver solution
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  if(!is_timedep(solvedExample))
    solve_exact(solvedExample, space, diffusivity, s, sigma, exact_solution, initial_sln, time_step_length, logger, logger_details);
  
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();
  if(algorithm == Multiscale || algorithm == Both)
  {
    logger.info("Multiscale solver");
    
    multiscale_decomposition(mesh, solvedExample, polynomialDegree, initialRefinementsCount, previous_mean_values, previous_derivatives, diffusivity, s, sigma, time_step_length,
    initial_sln, solution, exact_solution, &solution_view, &exact_view, logger, logger_details);
    
    cpu_time.tick();
    logger.info("Multiscale total: %f", cpu_time.last());
    logger.info("\n");
  }
  
  if(algorithm == pMultigrid || algorithm == Both)
  {
    Hermes::vector<int> steps(1, 3, 5, 10, 15);
    for(int si = 0; si < steps.size(); si++)
    {
      cpu_time.tick();
      logger.info("p-Multigrid solver - %i steps", steps[si]);

      MeshFunctionSharedPtr<double> previous_solution_local(new ExactSolutionMovingPeak(mesh, MovingPeakDiffusivity, M_PI / 2.));
  
      p_multigrid(mesh, solvedExample, polynomialDegree, initialRefinementsCount, previous_solution_local, diffusivity, time_step_length,
        solution, exact_solution, &solution_view, &exact_view, s, sigma, logger, logger_details, steps[si]);
        
      cpu_time.tick();
      logger.info("p-Multigrid solver - %i steps total: %f", steps[si], cpu_time.last());
    }
  }
  
  bool onlySmoothing = false;
  if(onlySmoothing)
  {
    cpu_time.tick();
    logger.info("only Smoothing solver");

    smoothing(mesh, solvedExample, polynomialDegree, previous_solution, diffusivity, time_step_length, 
      solution, exact_solution, &solution_view, &exact_view, s, sigma, logger, logger_details, 1);
      
    cpu_time.tick();
    logger.info("only Smoothing total: %s", cpu_time.last_str().c_str());
  }
  
  //View::wait();
  return 0;
}
