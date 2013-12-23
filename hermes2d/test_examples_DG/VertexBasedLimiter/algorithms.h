#include "hermes2d.h"
#include "definitions.h"

// Under relaxation in Multiscale
#define OMEGA 1.0

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;
using namespace Hermes::Algebra;

bool is_timedep(SolvedExample solvedExample);

void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step, int poly_degree, int init_ref_num);

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values,
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length, MeshFunctionSharedPtr<double> previous_solution,  MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, double cfl, int step);

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_sln,
                              double diffusivity, double time_step_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, int steps, double cfl, int V_cycles_per_time_step);
                              
void smoothing(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                              double diffusivity, double time_step_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, int steps);
