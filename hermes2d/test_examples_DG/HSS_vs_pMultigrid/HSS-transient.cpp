for (int step = 1; step <= iteration_count; step++)
{
  static_log.info("\tIteration: %i.", step);
  iterations++;

  // 1. predictor
  // M
  if (is_timedep(solvedExample))
    matrix_M_der.multiply_with_vector(sln_der.v, vector_A_der.v, true);
  else
    matrix_M_der.multiply_with_vector(sln_der_k.v, vector_A_der.v, true);
  add_ders(&sln_means_k, &sln_means_long, const_space, full_space);
  // -B
  matrix_A_full.multiply_with_vector(sln_means_long.v, sln_means_long_temp.v, true);
  SimpleVector<double>* temp_2 = (SimpleVector<double>*)cut_off_means(sln_means_long_temp.v, space, full_space)->change_sign();
  vector_A_der.add_vector(temp_2);
  delete temp_2;
  // (A-A~) - offdiag
  matrix_A_offdiag.multiply_with_vector(sln_der_k.v, sln_der_offdiag.v, true);
  vector_A_der.add_vector(sln_der_offdiag.change_sign());
  // b'
  vector_A_der.add_vector(&vector_b_der);
  // SOLVE
  solver_der.solve();
  sln_der_k_tilde.set_vector(solver_der.get_sln_vector());

  // 2. means
  // M
  if (is_timedep(solvedExample))
    matrix_M_means.multiply_with_vector(sln_means.v, vector_A_means.v, true);
  else
    matrix_M_means.multiply_with_vector(sln_means_k.v, vector_A_means.v, true);
  // -B
  add_means(&sln_der_k_tilde, &sln_der_long, space, full_space);
  matrix_A_full.multiply_with_vector(sln_der_long.v, sln_der_long_temp.v, true);
  SimpleVector<double>* temp_1 = (SimpleVector<double>*)cut_off_ders(sln_der_long_temp.v, const_space, full_space)->change_sign();
  vector_A_means.add_vector(temp_1);
  delete temp_1;
  // b
  vector_A_means.add_vector(&vector_b_means);
  // SOLVE
  solver_means.solve();
  sln_means_k.set_vector(solver_means.get_sln_vector());

  // 3. corrector
  // M
  if (is_timedep(solvedExample))
    matrix_M_der.multiply_with_vector(sln_der.v, vector_A_der.v, true);
  else
    matrix_M_der.multiply_with_vector(sln_der_k.v, vector_A_der.v, true);
  // -B
  add_ders(&sln_means_k, &sln_means_long, const_space, full_space);
  matrix_A_full.multiply_with_vector(sln_means_long.v, sln_means_long_temp.v, true);
  SimpleVector<double>* temp_3 = (SimpleVector<double>*)cut_off_means(sln_means_long_temp.v, space, full_space)->change_sign();
  vector_A_der.add_vector(temp_3);
  delete temp_3;
  // (A - A~) - offdiag
  matrix_A_offdiag.multiply_with_vector(sln_der_k.v, sln_der_offdiag.v, true);
  vector_A_der.add_vector(sln_der_offdiag.change_sign());
  // b
  vector_A_der.add_vector(&vector_b_der);
  // SOLVE
  solver_der.solve();
  if (OMEGA >= 0.99)
    sln_der_k.set_vector(solver_der.get_sln_vector());
  else
  {
    for (int i = 0; i < ndofs; i++)
      sln_der_k.set(i, (OMEGA * solver_der.get_sln_vector()[i]) + ((1. - OMEGA) * sln_der_k.get(i)));
  }
  merged_sln = merge_slns(sln_means_k.v, const_space, sln_der_k.v, space, full_space);
}