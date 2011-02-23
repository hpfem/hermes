//////////////////////////// SDIRK-22 ////////////
template<typename Real, typename Scalar>
Scalar jac_sdirk(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / time_step
                       + BUTCHER_A_11 * (dlam_du(u_prev_newton->val[i]) * u->val[i] 
                                  * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                                  + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])));
  return result;
}

template<typename Real, typename Scalar>
Scalar res_sdirk_stage_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y1_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (Y1_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / time_step;
  }
  result += GAMMA * res_ss(n, wt, u_ext, v, e, ext, current_time + BUTCHER_C_1 * time_step);
  return result;
}

template<typename Real, typename Scalar>
Scalar res_ss(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext, double t)
{
  Scalar result = 0;
  Func<Scalar>* Y_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (lam(Y_prev_newton->val[i]) * (Y_prev_newton->dx[i] * v->dx[i] + Y_prev_newton->dy[i] * v->dy[i])
                       - heat_src(e->x[i], e->y[i], t) * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar res_sdirk_stage_2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y2_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  Func<Scalar>* Y1[] = {ext->fn[1]};
  for (int i = 0; i < n; i++) {
    result += wt[i] * (Y2_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / time_step;
  }
  result += BUTCHER_B_1 * res_ss(n, wt, Y1, v, e, ext, current_time + BUTCHER_C_1 * time_step);
  result += BUTCHER_B_2 * res_ss(n, wt, u_ext, v, e, ext, current_time + BUTCHER_C_2 * time_step);
  return result;
}

//////////////////////////// SDIRK-22 OLD FASHION ////////////
// Y1 ///////////
template<typename Real, typename Scalar>
Scalar jac1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y1_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]
                       + GAMMA * time_step * (dlam_du(Y1_prev_newton->val[i]) * u->val[i] * (Y1_prev_newton->dx[i] * v->dx[i] 
                                                                                       + Y1_prev_newton->dy[i] * v->dy[i])
                                        + lam(Y1_prev_newton->val[i]) * (u->dx[i] * v->dx[i] 
                                                                         + u->dy[i] * v->dy[i])));
  return result;
}

template<typename Real, typename Scalar>
Scalar res1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y1_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((Y1_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] 
                       + GAMMA * time_step * (lam(Y1_prev_newton->val[i]) * (Y1_prev_newton->dx[i] * v->dx[i] 
                                                                       + Y1_prev_newton->dy[i] * v->dy[i])
                                        - heat_src(e->x[i], e->y[i], current_time + GAMMA*time_step) * v->val[i]));
  return result;
}

// Y2 ///////////
template<typename Real, typename Scalar>
Scalar jac2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y2_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]
                       + GAMMA * time_step * (dlam_du(Y2_prev_newton->val[i]) * u->val[i] * (Y2_prev_newton->dx[i] * v->dx[i]
                                                                                       + Y2_prev_newton->dy[i] * v->dy[i])
                                        + lam(Y2_prev_newton->val[i]) * (u->dx[i] * v->dx[i]
                                                                         + u->dy[i] * v->dy[i])));
  return result;
}

template<typename Real, typename Scalar>
Scalar res2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y2_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  Func<Scalar>* Y1 = ext->fn[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((Y2_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] 
                       - (1-GAMMA)/GAMMA * (Y1->val[i] - u_prev_time->val[i]) * v->val[i]
                       + GAMMA * time_step * (lam(Y2_prev_newton->val[i]) * (Y2_prev_newton->dx[i] * v->dx[i]
                                                                 + Y2_prev_newton->dy[i] * v->dy[i])
                                        - heat_src(e->x[i], e->y[i], current_time + time_step) * v->val[i]));
  return result;
}

