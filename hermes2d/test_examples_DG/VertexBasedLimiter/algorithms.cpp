#include "algorithms.h"
static std::string SolvedExampleString[5] = { "1D", "CircularConvection", "MovingPeak", "AdvectedCube", "SolidBodyRotation" };
static double exact_solver_error;
static const double tolerance = 1e-4;
double initial_error = -1;
MeshFunctionSharedPtr<double> es(new Solution<double>());
double* es_v;

static Hermes::Mixins::Loggable static_log(true);

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

double calc_l2_error(SolvedExample solvedExample, MeshSharedPtr mesh, MeshFunctionSharedPtr<double> fn_1, MeshFunctionSharedPtr<double> fn_2, Hermes::Mixins::Loggable& logger)
{
  ErrorWeakForm wf(solvedExample);
  SpaceSharedPtr<double> mspace(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  wf.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(fn_1, fn_2));
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, mspace);
  SimpleVector<double> vector;
  dp->assemble(&vector);
  double result = 0.;
  for (int i = 0; i < vector.get_size(); i++)
    result += vector.get(i);
  result = std::sqrt(result);
  //logger.info("L2 Error: %g.", result);
  delete dp;
  return result;
}

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

bool error_condition(double error)
{
  return std::abs(error) < tolerance;
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

  ss_bmp << "solution_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".bmp";
  ss_vtk << "solution_" << SolvedExampleString[solvedExample] << "_" << init_ref_num << "_" << diffusivity << ".dat";
  // exact_solver_view->show(es);
  exact_solver_view->save_screenshot(ss_bmp.str().c_str(), true);
  exact_solver_view->get_linearizer()->save_solution_tecplot(es, ss_vtk.str().c_str(), "solution");
}

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_mean_values,
  MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length,
  MeshFunctionSharedPtr<double> previous_solution, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
  ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, double cfl, int step)
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
  SimpleVector<double> sln_der_long(full_ndofs);
  SimpleVector<double> sln_der_long_temp(full_ndofs);
  SimpleVector<double> sln_der_offdiag(ndofs);

  OGProjection<double>::project_global(const_space, previous_mean_values, sln_means.v);
  OGProjection<double>::project_global(space, previous_derivatives, sln_der.v);

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;

  double* merged_sln;

  double time = 0.;
  int iteration_count = (int)(is_timedep(solvedExample) ? std::ceil(end_time(solvedExample) / time_step_length) : 10000);
  for (int iteration = 0; iteration < iteration_count; iteration++)
  {
    iterations++;
    num_coarse++;
    static_log.info("Time step: %i, time: %f.", iteration, time);

    for (int iter_i = 0; iter_i < step; iter_i++)
    {
      matrix_M_means.multiply_with_vector(sln_means.v, vector_A_means.v, true);
      if (polynomialDegree)
      {
        add_means(&sln_der_k, &sln_der_long, space, full_space);
        matrix_A_full.multiply_with_vector(sln_der_long.v, sln_der_long_temp.v, true);
        SimpleVector<double>* temp_1 = (SimpleVector<double>*)cut_off_ders(sln_der_long_temp.v, const_space, full_space)->change_sign();
        vector_A_means.add_vector(temp_1);
        delete temp_1;
      }
      vector_A_means.add_vector(&vector_b_means);
      solver_means.solve();
      sln_means_k.set_vector(solver_means.get_sln_vector());

      if (polynomialDegree)
      {
        num_fine++;

        matrix_M_der.multiply_with_vector(sln_der.v, vector_A_der.v, true);
        add_ders(&sln_means_k, &sln_means_long, const_space, full_space);
        matrix_A_full.multiply_with_vector(sln_means_long.v, sln_means_long_temp.v, true);
        SimpleVector<double>* temp_2 = (SimpleVector<double>*)cut_off_means(sln_means_long_temp.v, space, full_space)->change_sign();
        vector_A_der.add_vector(temp_2);
        delete temp_2;

        matrix_A_offdiag.multiply_with_vector(sln_der_k.v, sln_der_offdiag.v, true);
        vector_A_der.add_vector(sln_der_offdiag.change_sign());

        vector_A_der.add_vector(&vector_b_der);
        solver_der.solve();
        if (iteration == 1 || OMEGA >= 0.99)
          sln_der_k.set_vector(solver_der.get_sln_vector());
        else
        {
          for (int i = 0; i < ndofs; i++)
            sln_der_k.set(i, (OMEGA * solver_der.get_sln_vector()[i]) + ((1. - OMEGA) * sln_der_k.get(i)));
        }
        merged_sln = merge_slns(sln_means_k.v, const_space, sln_der_k.v, space, full_space);
      }
      else
      {
        merged_sln = sln_means.v;
      }
    }
    sln_means.set_vector(&sln_means_k);
    sln_der.set_vector(&sln_der_k);
    /*
    if (polynomialDegree)
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    else
      Solution<double>::vector_to_solution(sln_means.v, const_space, solution);

    solution_view->show(solution);

    if (iteration == iteration_count - 1)
    {
    
    std::stringstream ss_bmp, ss_vtk;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_p=" << polynomialDegree << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << ".dat";
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution");
    }
    */

    bool done = !is_timedep(solvedExample) && error_reduction_condition(calc_l2_error_algebraic(polynomialDegree ? full_space : const_space, merged_sln, es_v, &logger_details, iteration, init_ref_num, diffusivity));

    if (is_timedep(solvedExample))
    {
      if (time + time_step_length > end_time(solvedExample))
      {
        time_step_length = end_time(solvedExample) - time;
        time = end_time(solvedExample);
      }
    }

    time += time_step_length;

    bool finish_timedep = is_timedep(solvedExample) && (iteration == iteration_count - 1);
    
    if (polynomialDegree && !finish_timedep)
      delete[] merged_sln;

    if (done)
      break;
  }

  if (is_timedep(solvedExample))
  {
    ((ExactSolutionMovingPeak*)exact_solution.get())->set_current_time(time + (M_PI / 2.));
    exact_view->show(exact_solution);
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

    if (polynomialDegree)
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
    else
      Solution<double>::vector_to_solution(sln_means.v, const_space, solution);

    std::stringstream ss_vtk;
    std::stringstream ss_bmp;
    std::stringstream ss_bmpe;
    std::stringstream ss_vtke;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_" << "HSS(" << step << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";
    ss_vtke.precision(2);
    ss_vtke.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtke << "exact_solution_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";

    ss_bmp.precision(2);
    ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmp << "solution_" << "HSS(" << step << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";
    ss_bmpe.precision(2);
    ss_bmpe.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmpe << "exact_solution_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";
    
    exact_view->show(exact_solution);
    //exact_view->save_screenshot(ss_bmpe.str().c_str(), true);
    
    solution_view->show(solution);
    //solution_view->save_screenshot(ss_bmp.str().c_str(), true);
    
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);
    solution_view->get_linearizer()->save_solution_tecplot(exact_solution, ss_vtke.str().c_str(), "exactSolution", 1, 2.0);

    errorCalculator.calculate_errors(solution, exact_solution);
    logger.info("%f", std::sqrt(errorCalculator.get_total_error_squared()));
  }
  logger.info("%i", iterations);
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, int init_ref_num, MeshFunctionSharedPtr<double> previous_sln,
  double diffusivity, double time_step_length,
  MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
  ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, int steps, double cfl, int V_cycles_per_time_step)
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

  double initial_residual_norm;
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
  int iteration_count = (int)(is_timedep(solvedExample) ? std::ceil(end_time(solvedExample) / time_step_length) : 1);
  for (int time_step = 0; time_step < iteration_count; time_step++)
  {
    if(!is_timedep(solvedExample))
      static_log.info("Time step: %i, time: %f.", time_step, time);

    for (int step = 0; step < (is_timedep(solvedExample) ? V_cycles_per_time_step : iteration_count); step++)
    {
      static_log.info("V-cycle %i.", step);
      v_cycles++;

#pragma region 0 - highest level
      // Store the previous solution.
      if (polynomialDegree > 1)
      {
        for (int iteration = 1; iteration <= steps; iteration++)
        {
          // Solve for increment.
          matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
          vector_A_2.change_sign()->add_vector(&vector_b_2);

          util_2.set_vector(&prev_sln_2);
          util_2.change_sign()->add_vector(sln_2.v)->change_sign();
          matrix_M_2.multiply_with_vector(util_2.v, util_21.v, true);

          vector_A_2.add_vector(util_21.v);

          solver_2.solve();
          sln_2.add_vector(solver_2.get_sln_vector());
        }
      }

#pragma endregion

#pragma region 1 - intermediate level

      // f_P1
      SimpleVector<double> f_P1(ndofs_1);
      f_P1.zero();
      // Minus A_P1
      SimpleVector<double> R_P1(ndofs_1);
      // Minus(minus) projected_A_P1
      SimpleVector<double> projected_A_2(ndofs_2);
      matrix_A_2.multiply_with_vector(sln_2.v, projected_A_2.v, true);

      SimpleVector<double>* projected_A_P_1
        = (SimpleVector<double>*)cut_off_quadratic_part(projected_A_2.v, space_1, space_2);

      SimpleVector<double>* sln_2_projected = cut_off_quadratic_part(sln_2.v, space_1, space_2);
      matrix_A_1.multiply_with_vector(sln_2_projected->v, R_P1.v, true);

      sln_1.set_vector(sln_2_projected);

      R_P1.change_sign();
      f_P1.add_vector(&R_P1);
      f_P1.add_vector(projected_A_P_1);
      delete projected_A_P_1;
      delete sln_2_projected;
      f_P1.change_sign();

      for (int iteration = 1; iteration <= steps; iteration++)
      {
        // Solve for increment.
        SimpleVector<double>* rhs;

        // A(u_K) - done after the first step.
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if (polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);

        util_1.set_vector(&prev_sln_1);
        util_1.change_sign()->add_vector(sln_1.v)->change_sign();
        matrix_M_1.multiply_with_vector(util_1.v, util_11.v, true);

        vector_A_1.add_vector(util_11.v);

        solver_1.solve();
        sln_1.add_vector(solver_1.get_sln_vector());
      }

#pragma endregion

#pragma region  2 - Solve the problem on the coarse level exactly

      // f_P0
      SimpleVector<double> f_P0(ndofs_0);
      f_P0.zero();
      // Minus A_P0
      SimpleVector<double> R_P0(ndofs_0);
      // Minus(minus) projected_A_P0
      SimpleVector<double> projected_A_1(ndofs_1);
      matrix_A_1.multiply_with_vector(sln_1.v, projected_A_1.v, true);

      SimpleVector<double>* projected_A_P_0
        = (SimpleVector<double>*)cut_off_linear_part(projected_A_1.v, space_0, space_1);

      SimpleVector<double>* sln_1_projected = cut_off_linear_part(sln_1.v, space_0, space_1);
      matrix_A_0.multiply_with_vector(sln_1_projected->v, R_P0.v, true);

      SimpleVector<double> projected_f_P1(ndofs_1);
      projected_f_P1.set_vector(&f_P1);

      sln_0.set_vector(sln_1_projected);

      R_P0.change_sign();
      f_P0.add_vector(&R_P0);
      f_P0.add_vector(projected_A_P_0);
      if (polynomialDegree > 1)
      {
        Vector<double>* temp = cut_off_linear_part(projected_f_P1.v, space_0, space_1)->change_sign();
        f_P0.add_vector(temp);
        delete temp;
      }
      delete projected_A_P_0;
      delete sln_1_projected;
      f_P0.change_sign();

      num_coarse++;

      // A(u_K) - done after the first step.
      matrix_M_0.multiply_with_vector(prev_sln_0.v, vector_A_0.v, true);
      vector_A_0.add_vector(&f_P0)->add_vector(&vector_b_0);

      solver_0.solve();
      sln_0.set_vector(solver_0.get_sln_vector());

#pragma endregion

#pragma region 1 - intermediate level
      // Store the previous solution.
      double* vector_ = merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, false);
      sln_1.set_vector(vector_);
      delete[] vector_;

      for (int iteration = 1; iteration <= steps; iteration++)
      {
        // Solve for increment.
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if (polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);

        util_1.set_vector(&prev_sln_1);
        util_1.change_sign()->add_vector(sln_1.v)->change_sign();
        matrix_M_1.multiply_with_vector(util_1.v, util_11.v, true);

        vector_A_1.add_vector(util_11.v);

        solver_1.solve();
        sln_1.add_vector(solver_1.get_sln_vector());
      }

#pragma endregion

#pragma region 0 - highest level

      if (polynomialDegree > 1)
      {
        // Store the previous solution.
        double* vector_ = merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false);
        sln_2.set_vector(vector_);
        delete[] vector_;

        for (int iteration = 1; iteration <= steps; iteration++)
        {
          // Solve for increment.
          matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
          vector_A_2.change_sign()->add_vector(&vector_b_2);
          
          util_2.set_vector(&prev_sln_2);
          util_2.change_sign()->add_vector(sln_2.v)->change_sign();
          matrix_M_2.multiply_with_vector(util_2.v, util_21.v, true);

          vector_A_2.add_vector(util_21.v);

          solver_2.solve();
          sln_2.add_vector(solver_2.get_sln_vector());
        }


        /*
        // Make solution & display.
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->set_title("Time: %f", time);
        solution_view->show(previous_sln);

        if ((step == 1) || step > iteration_count - 2)
        {
        std::stringstream ss_bmp, ss_vtk;
        ss_bmp.precision(2);
        ss_vtk.precision(2);
        ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
        ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);

        ss_bmp << "solution_p=" << polynomialDegree << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << ".bmp";
        ss_vtk << "solution_p=" << polynomialDegree << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << ".dat";
        solution_view->save_screenshot(ss_bmp.str().c_str(), true);
        solution_view->get_linearizer()->save_solution_tecplot(previous_sln, ss_vtk.str().c_str(), "solution");
        }
        */
      }
#pragma endregion

      if (!is_timedep(solvedExample) && error_reduction_condition(calc_l2_error_algebraic(space_2, sln_2.v, es_v)))
        break;
    }

    if (is_timedep(solvedExample))
    {
      if (time + time_step_length > end_time(solvedExample))
      {
        time_step_length = end_time(solvedExample) - time;
        time = end_time(solvedExample);
      }
    }

    time += time_step_length;

    prev_sln_2.set_vector(&sln_2);
    prev_sln_1.set_vector(&sln_1);
    prev_sln_0.set_vector(&sln_0);
  }

  if (is_timedep(solvedExample))
  {
    ((ExactSolutionMovingPeak*)exact_solution.get())->set_current_time(time + (M_PI / 2.));
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);

    Solution<double>::vector_to_solution(&sln_2, space_2, solution);
    
    std::stringstream ss_vtk;
    std::stringstream ss_vtke;
    ss_vtk.precision(2);
    ss_vtk.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_vtk << "solution_" << "MG(" << V_cycles_per_time_step << "-" << steps << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".dat";
    solution_view->get_linearizer()->save_solution_tecplot(solution, ss_vtk.str().c_str(), "solution", 1, 2.0);

    std::stringstream ss_bmp;
    std::stringstream ss_bmpe;
    ss_bmp.precision(2);
    ss_bmp.setf(std::ios_base::uppercase | std::ios_base::scientific);
    ss_bmp << "solution_" << "MG(" << V_cycles_per_time_step << "-" << steps << ")_" << SolvedExampleString[solvedExample] << "_meshRefs=" << init_ref_num << "_D=" << diffusivity << "_CFL=" << cfl << ".bmp";
    
    solution_view->show(solution);
    solution_view->save_screenshot(ss_bmp.str().c_str(), true);

    errorCalculator.calculate_errors(solution, exact_solution);
    logger.info("%f", std::sqrt(errorCalculator.get_total_error_squared()));
  }
  logger.info("%i", v_cycles);
}

void smoothing(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
  double diffusivity, double time_step_length,
  MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
  ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, Hermes::Mixins::Loggable& logger_details, int steps)
{
  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_2;
  SimpleVector<double> vector_b_2;

  // Matrices (M+A_tilde), vectors -A(u_K)
  SmoothingWeakForm weakform_smoother(solvedExample, false, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  SmoothingWeakForm weakform_smoother_coarse(solvedExample, false, 1, add_inlet(solvedExample), "Inlet", diffusivity, s, sigma);
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  weakform_smoother_coarse.set_current_time_step(time_step_length);
  weakform_smoother_coarse.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2(ndofs_2);

  // Assembler.
  DiscreteProblem<double> dp;
  // Level 2.
  dp.set_space(space_2);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_2, &vector_b_2);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_2);

  UMFPackLinearMatrixSolver<double> solver_2(&matrix_MA_tilde_2, &vector_A_2);
  solver_2.setup_factorization();
  solver_2.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  double* residual_2 = new double[ndofs_2];
  SimpleVector<double> sln_2(ndofs_2);
  sln_2.zero();

  double initial_residual_norm;

  // Reports.
  int v_cycles = 0;

  for (int step = 1;; step++)
  {
    logger_details.info("V-cycle %i.", step);
    v_cycles++;

    // Solve for increment.
    matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
    vector_A_2.change_sign()->add_vector(&vector_b_2);
    solver_2.solve();
    sln_2.add_vector(solver_2.get_sln_vector());

    // Make solution
    Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);

    // Residual check.
    // residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger_details, step, true);

    // Error & exact solution display.
    //solution_view->show(->show(previous_sln);

    if (error_condition(calc_l2_error(solvedExample, mesh, previous_sln, es, logger_details)))
      break;
  }
  logger.info("%i", v_cycles);
}
