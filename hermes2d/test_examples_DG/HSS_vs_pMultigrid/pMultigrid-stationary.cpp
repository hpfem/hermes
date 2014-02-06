for (int step = 0; step < iteration_count; step++)
{
  static_log.info("V-cycle %i.", step);
  v_cycles++;
  static_log.info("Time step: %i, time: %f.", step, time);

#pragma region 0 - highest level
  // Store the previous solution.
  if (polynomialDegree > 1)
  {
    for (int smoothing_step = 1; smoothing_step <= smoothing_steps_per_V_cycle; smoothing_step++)
    {
      // Solve for increment.
      matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
      vector_A_2.change_sign()->add_vector(&vector_b_2);
      solver_2.solve();
      sln_2.add_vector(solver_2.get_sln_vector());
    }

    // Residual calculation.
    matrix_A_2.multiply_with_vector(sln_2.v, residual_2, true);
    for (int i = 0; i < matrix_A_2.get_size(); i++)
      residual_2[i] -= vector_b_2.get(i);
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
    if (iteration == 1)
    {
      if (polynomialDegree > 1)
      {
        rhs = cut_off_quadratic_part(residual_2, space_1, space_2);
        memcpy(vector_A_1.v, rhs->v, ndofs_1 * sizeof(double));
        delete rhs;
      }
      else
      {
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);
        vector_A_1.change_sign()->add_vector(&vector_b_1);
      }
    }
    else
    {
      // A(u_K) - done after the first step.
      matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

      if (polynomialDegree > 1)
        vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
      else
        vector_A_1.change_sign()->add_vector(&vector_b_1);
    }
    solver_1.solve();
    sln_1.add_vector(solver_1.get_sln_vector());
  }

  // Residual calculation.
  matrix_A_1.multiply_with_vector(sln_1.v, residual_1, true);
  for (int i = 0; i < matrix_A_1.get_size(); i++)
    residual_1[i] -= vector_b_1.get(i);
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
  matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);
  vector_A_0.change_sign()->add_vector(&f_P0)->add_vector(&vector_b_0);

  solver_0.solve();
  sln_0.add_vector(solver_0.get_sln_vector());
#pragma endregion

#pragma region 1 - intermediate level
  // Store the previous solution.
  sln_1.set_vector(merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, false));

  for (int iteration = 1; iteration <= steps; iteration++)
  {
    // Solve for increment.
    matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

    if (polynomialDegree > 1)
      vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
    else
      vector_A_1.change_sign()->add_vector(&vector_b_1);

    solver_1.solve();
    sln_1.add_vector(solver_1.get_sln_vector());
  }

  // Residual calculation.
  matrix_A_1.multiply_with_vector(sln_1.v, residual_1, true);
  for (int i = 0; i < matrix_A_1.get_size(); i++)
    residual_1[i] -= vector_b_1.get(i);
#pragma endregion

#pragma region 0 - highest level

  if (polynomialDegree > 1)
  {
    // Store the previous solution.
    sln_2.set_vector(merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false));

    for (int iteration = 1; iteration <= steps; iteration++)
    {
      // Solve for increment.
      matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
      vector_A_2.change_sign()->add_vector(&vector_b_2);
      solver_2.solve();
      sln_2.add_vector(solver_2.get_sln_vector());
    }

    // Residual calculation.
    matrix_A_2.multiply_with_vector(sln_2.v, residual_2, true);
    for (int i = 0; i < matrix_A_2.get_size(); i++)
      residual_2[i] -= vector_b_2.get(i);


    /*
    if ((step == 1) || step > iteration_count - 2)
    {
    // Make solution & display.
    Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
    solution_view->set_title("Time: %f", time);
    solution_view->show(previous_sln);

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

  if (is_timedep(solvedExample))
  {
    if (time + time_step_length > end_time(solvedExample))
    {
      time_step_length = end_time(solvedExample) - time;
      time = end_time(solvedExample);
    }
  }

  time += time_step_length;

  if (!is_timedep(solvedExample) && error_reduction_condition(calc_l2_error_algebraic(space_2, sln_2.v, es_v)))
    break;
}