#include "algorithms.h"

// Relative tolerance drop (1e-4 == 4 orders of magnitude drop)
static const double tolerance = 1e-4;

static const int integrationOrder = 4;

// Utilities
static std::string SolvedExampleString[5] = { "1D", "CircularConvection", "MovingPeak", "AdvectedCube", "SolidBodyRotation" };
static double exact_solver_error;
double initial_error;
MeshFunctionSharedPtr<double> es(new Solution<double>());
double* es_v;

// Uncomment to have OpenGL output throughout calculation.
#define SHOW_OUTPUT

// Under relaxation in Multiscale
#define OMEGA 1.0

// Static logging for output in terminal.
static Hermes::Mixins::Loggable static_log(true);

double calc_l2_error_algebraic(SpaceSharedPtr<double> space, double* v1, double* v2, Hermes::Mixins::Loggable* logger = NULL, int iteration = 0, int init_refs = 0, double D = 0.)
{
  double result = 0.;
  for (int i = 0; i < space->get_num_dofs(); i++)
    result += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  result = std::sqrt(result);
  if (logger)
    logger->info("%d,%f,%d,%f", init_refs, D, iteration, result);
  return result;
}

bool error_reduction_condition(double error)
{
  return std::abs(error / initial_error) < tolerance;
}

void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step, int poly_degree, int init_ref_num)
{
  MeshFunctionSharedPtr<double> exact_solver_sln(new Solution<double>());
  ScalarView* exact_solver_view = new ScalarView("Exact solver solution", new WinGeom(0, 360, 600, 350));

  // Exact solver
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, initial_sln);
  weakform_exact.set_current_time_step(time_step);
  LinearSolver<double> solver_exact(&weakform_exact, space);

  // Solve
  solver_exact.solve();
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, exact_solver_sln);

  // Initial error
  initial_error = get_l2_norm(solver_exact.get_sln_vector(), space->get_num_dofs());

  // es and es_v
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, es);
  es_v = new double[space->get_num_dofs()];
  memcpy(es_v, solver_exact.get_sln_vector(), sizeof(double)* space->get_num_dofs());

  // output
  std::stringstream ss_bmp, ss_vtk;
  ss_bmp.precision(2);
  ss_vtk.precision(2);
  ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);

  ss_bmp << "exact_solution_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".bmp";
  ss_vtk << "exact_solution_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".dat";
#ifdef SHOW_OUTPUT
  exact_solver_view->show(es);
  exact_solver_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
  exact_solver_view->get_linearizer()->save_solution_tecplot(es, ss_vtk.str().c_str(), "solution");
}

std::string multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values,
  MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length,
  MeshFunctionSharedPtr<double> previous_solution, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
  ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, double cfl, int steps_per_time_step)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));

  int ndofs = space->get_num_dofs();
  int const_ndofs = const_space->get_num_dofs();
  int full_ndofs = full_space->get_num_dofs();

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  MultiscaleWeakForm weakform_implicit(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution, false);
  MultiscaleWeakForm weakform_explicit(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution, true);
  ExplicitWeakFormOffDiag weakform_explicit_offdiag(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  MassWeakForm weakform_mass;
  weakform_exact.set_current_time_step(time_step_length);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_explicit.set_current_time_step(time_step_length);
  weakform_explicit_offdiag.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_full;
  CSCMatrix<double> matrix_A_means_just_A;
  CSCMatrix<double> matrix_A_der;
  SimpleVector<double> vector_b_der;
  CSCMatrix<double> matrix_M_der;
  CSCMatrix<double> matrix_A_offdiag;

  CSCMatrix<double> matrix_A_means;
  SimpleVector<double> vector_b_means;
  CSCMatrix<double> matrix_M_means;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_full);
  dp.set_space(const_space);
  dp.assemble(&matrix_A_means_just_A);

  // Level 1.
  dp.set_space(space);
  dp.set_weak_formulation(&weakform_explicit);
  dp.assemble(&matrix_A_der, &vector_b_der);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_der);
  dp.set_weak_formulation(&weakform_explicit_offdiag);
  dp.assemble(&matrix_A_offdiag);

  // Level 0.
  dp.set_space(const_space);
  dp.set_weak_formulation(&weakform_implicit);
  dp.assemble(&matrix_A_means, &vector_b_means);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_means);

  SimpleVector<double> vector_A_der(ndofs);
  SimpleVector<double> vector_A_means(const_ndofs);

  UMFPackLinearMatrixSolver<double> solver_means(&matrix_A_means, &vector_A_means);
  solver_means.setup_factorization();
  solver_means.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_der(&matrix_A_der, &vector_A_der);
  solver_der.setup_factorization();
  solver_der.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  SimpleVector<double> sln_means(const_ndofs);
  SimpleVector<double> sln_means_k(const_ndofs);
  SimpleVector<double> sln_means_long(full_ndofs);
  SimpleVector<double> sln_means_long_temp(full_ndofs);
  SimpleVector<double> sln_der(ndofs);
  SimpleVector<double> sln_der_k(ndofs);
  SimpleVector<double> sln_der_k_tilde(ndofs);
  SimpleVector<double> sln_der_long(full_ndofs);
  SimpleVector<double> sln_der_long_temp(full_ndofs);
  SimpleVector<double> sln_der_offdiag(ndofs);

  OGProjection<double>::project_global(const_space, previous_mean_values, sln_means.v);
  OGProjection<double>::project_global(space, previous_derivatives, sln_der.v);

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;
  double time = 0.;

  double* merged_sln;

  int time_step_count = (int)(is_timedep(solvedExample) ? std::ceil(end_time(solvedExample) / time_step_length) : 10000);
  int iteration_count = steps_per_time_step;
  for (int time_step = 1; time_step <= time_step_count; time_step++)
  {
    if(is_timedep(solvedExample))
      static_log.info("Time step: %i, time: %f.", time_step, time);
    else
      static_log.info("Time step: %i.", time_step);

    if (is_timedep(solvedExample))
    {
#include "HSS-transient.cpp"
    }
    else
    {
#include "HSS-stationary.cpp"
    }
    sln_means.set_vector(&sln_means_k);
    sln_der.set_vector(&sln_der_k);
#ifdef SHOW_OUTPUT
    if (polynomialDegree)
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    else
      Solution<double>::vector_to_solution(sln_means.v, const_space, solution);

    solution_view->set_title("Time: %f.");
    solution_view->show(solution);
#endif

    bool done = !is_timedep(solvedExample) && error_reduction_condition(calc_l2_error_algebraic(polynomialDegree ? full_space : const_space, merged_sln, es_v, &logger_details, time_step, init_ref_num, diffusivity));

    if (is_timedep(solvedExample))
    {
      if (time + time_step_length > end_time(solvedExample))
      {
        time_step_length = end_time(solvedExample) - time;
        time = end_time(solvedExample);
      }
      else
        time += time_step_length;
    }

    bool finish_timedep = is_timedep(solvedExample) && (time_step == time_step_count);

    if (polynomialDegree && !finish_timedep)
      delete[] merged_sln;

    if (done)
      break;
  }

  std::stringstream outStream;
  outStream << iterations;
  if (is_timedep(solvedExample))
  {
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

    if (polynomialDegree)
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    else
      Solution<double>::vector_to_solution(sln_means.v, const_space, solution);

    std::stringstream ss_vtk;
    std::stringstream ss_bmp;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_" << "HSS(" << steps_per_time_step << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";

    ss_bmp.precision(2);
    ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmp << "solution_" << "HSS(" << steps_per_time_step << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
    solution_view->show(solution);
    solution_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);

    errorCalculator.calculate_errors(solution, es);

    outStream << "|" << std::sqrt(errorCalculator.get_total_error_squared());
  }

  return outStream.str();
}

std::string p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_sln,
  double diffusivity, double time_step_length,
  MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
  ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, int smoothing_steps_per_V_cycle, double cfl, int V_cycles_per_time_step)
{
  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_2;
  SimpleVector<double> vector_b_2;
  CSCMatrix<double> matrix_A_1;
  SimpleVector<double> vector_b_1;
  CSCMatrix<double> matrix_A_0;
  SimpleVector<double> vector_b_0;

  // Matrices (M+A_tilde), vectors -A(u_K)
  SmoothingWeakForm weakform_smoother(solvedExample, true, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  SmoothingWeakForm weakform_smoother_coarse(solvedExample, false, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  weakform_smoother_coarse.set_current_time_step(time_step_length);
  weakform_smoother_coarse.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  MassWeakForm weakform_mass;

  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2(ndofs_2);
  CSCMatrix<double> matrix_MA_tilde_1;
  SimpleVector<double> vector_A_1(ndofs_1);
  CSCMatrix<double> matrix_MA_0;
  SimpleVector<double> vector_A_0(ndofs_0);

  CSCMatrix<double> matrix_M_2;
  CSCMatrix<double> matrix_M_1;
  CSCMatrix<double> matrix_M_0;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(space_2);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_2, &vector_b_2);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_2);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_2);

  // Level 1.
  dp.set_space(space_1);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_1, &vector_b_1);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_1);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_1);

  // Level 0.
  dp.set_space(space_0);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_0, &vector_b_0);
  dp.set_weak_formulation(&weakform_smoother_coarse);
  dp.assemble(&matrix_MA_0);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_0);

  UMFPackLinearMatrixSolver<double> solver_2(&matrix_MA_tilde_2, &vector_A_2);
  solver_2.setup_factorization();
  solver_2.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_1(&matrix_MA_tilde_1, &vector_A_1);
  solver_1.setup_factorization();
  solver_1.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_0(&matrix_MA_0, &vector_A_0);
  solver_0.setup_factorization();
  solver_0.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  SimpleVector<double> sln_2(ndofs_2);
  SimpleVector<double> sln_1(ndofs_1);
  SimpleVector<double> sln_0(ndofs_0);

  SimpleVector<double> prev_sln_2(ndofs_2);
  SimpleVector<double> prev_sln_1(ndofs_1);
  SimpleVector<double> prev_sln_0(ndofs_0);

  SimpleVector<double> util_2(ndofs_2), util_21(ndofs_2);
  SimpleVector<double> util_1(ndofs_1), util_11(ndofs_2);
  SimpleVector<double> util_0(ndofs_0), util_01(ndofs_2);

  // Reports.
  int num_coarse = 0;
  int num_2 = 0;
  int num_1 = 0;
  int v_cycles = 0;

  OGProjection<double>::project_global(space_2, previous_sln, &sln_2);
  prev_sln_2.set_vector(&sln_2);
  prev_sln_1.set_vector(cut_off_quadratic_part(prev_sln_2.v, space_1, space_2));
  prev_sln_0.set_vector(cut_off_linear_part(prev_sln_1.v, space_0, space_1));

  double time = 0.;
  int time_step_count = (int)(is_timedep(solvedExample) ? std::ceil(end_time(solvedExample) / time_step_length) : 1);
  int iteration_count = (int)(is_timedep(solvedExample) ? V_cycles_per_time_step : 10000);
  for (int time_step = 1; time_step <= time_step_count; time_step++)
  {
    if(is_timedep(solvedExample))
      static_log.info("Time step: %i, time: %f.", time_step, time);
    else
      static_log.info("Time step: %i.", time_step);

    if (is_timedep(solvedExample))
    {
#include "pMultigrid-transient.cpp"
    }
    else
    {
#include "pMultigrid-stationary.cpp"
    }
    

    if (is_timedep(solvedExample))
    {
      if (time + time_step_length > end_time(solvedExample))
      {
        time_step_length = end_time(solvedExample) - time;
        time = end_time(solvedExample);
      }
      else
        time += time_step_length;
    }

    prev_sln_2.set_vector(&sln_2);
    prev_sln_1.set_vector(&sln_1);
    prev_sln_0.set_vector(&sln_0);
  }

  std::stringstream outStream;
  outStream << v_cycles;
  if (is_timedep(solvedExample))
  {
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

    Solution<double>::vector_to_solution(&sln_2, space_2, solution);

    std::stringstream ss_vtk;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_" << "MG(" << V_cycles_per_time_step << "-" << smoothing_steps_per_V_cycle << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);

    std::stringstream ss_bmp;
    ss_bmp.precision(2);
    ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmp << "solution_" << "MG(" << V_cycles_per_time_step << "-" << smoothing_steps_per_V_cycle << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
    solution_view->show(solution);
    solution_view->save_screenshot(ss_bmp.str().c_str(), true);
#endif
    errorCalculator.calculate_errors(solution, es);
    outStream << "|" << std::sqrt(errorCalculator.get_total_error_squared());
  }

  return outStream.str();
}

void exact_solver_timedep(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, double diffusivity, double s, double sigma, double time_step_length,
  MeshFunctionSharedPtr<double> previous_solution, MeshFunctionSharedPtr<double> exact_solution, ScalarView* exact_view, double cfl)
{
  if (!is_timedep(solvedExample))
    return;

  // Standard L2 space.
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int full_ndofs = full_space->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakFormTimedep weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  weakform_exact.set_ext(previous_solution);
  CSCMatrix<double> matrix_A;
  SimpleVector<double> vector_b;

  // Assembler.
  DiscreteProblem<double> dp;
  dp.set_global_integration_order(integrationOrder);
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A);

  UMFPackLinearMatrixSolver<double> solver(&matrix_A, &vector_b);
  solver.setup_factorization();
  solver.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Reporting.

  double time = 0.;
  int iteration_count = std::ceil(end_time(solvedExample) / time_step_length);
  for (int iteration = 0; iteration <= iteration_count; ++iteration)
  {
    dp.assemble(&vector_b);
    solver.solve();
    Solution<double>::vector_to_solution(solver.get_sln_vector(), full_space, previous_solution);
    static_log.info("Time step: %i, time: %f.", iteration, time);

#ifdef SHOW_OUTPUT
    exact_view->show(previous_solution);
#endif

    if (is_timedep(solvedExample))
    {
      if (time + time_step_length > end_time(solvedExample))
      {
        time_step_length = end_time(solvedExample) - time;
        time = end_time(solvedExample);
      }
      else
        time += time_step_length;
    }
  }
  Solution<double>::vector_to_solution(solver.get_sln_vector(), full_space, es);

  std::stringstream ss_bmpe;
  std::stringstream ss_vtke;
  ss_vtke.precision(2);
  ss_vtke.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_vtke << "exact_solution_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";

  ss_bmpe.precision(2);
  ss_bmpe.setf(std::ios_base::uppercase | std::ios_base::scientific);
  ss_bmpe << "exact_solution_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";

#ifdef SHOW_OUTPUT
  exact_view->show(es);
  exact_view->save_screenshot(ss_bmpe.str().c_str(), true);
#endif
  exact_view->get_linearizer()->save_solution_tecplot(es, ss_vtke.str().c_str(), "exactSolution", 1, 2.0);
}

// Utilities.
bool add_inlet(SolvedExample solvedExample)
{
  switch (solvedExample)
  {
  case SolidBodyRotation:
  case AdvectedCube:
  case MovingPeak:
    return false;
  case CircularConvection:
  case Benchmark:
    return true;
  }
}

bool is_timedep(SolvedExample solvedExample)
{
  switch (solvedExample)
  {
  case CircularConvection:
  case Benchmark:
    return false;
  case MovingPeak:
  case AdvectedCube:
  case SolidBodyRotation:
    return true;
  }
}

double end_time(SolvedExample solvedExample)
{
  switch (solvedExample)
  {
  case CircularConvection:
  case Benchmark:
    return 9999999999.;
  case MovingPeak:
    return M_PI * 2.;
  case AdvectedCube:
    return 1.;
  case SolidBodyRotation:
    return M_PI * 2.;
  }
}
