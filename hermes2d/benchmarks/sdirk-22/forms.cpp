// Jacobian matrix.
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] 
                       + TAU * (dlam_du(u_prev_newton->val[i]) * u->val[i] * (u_prev_newton->dx[i] * v->dx[i] 
                                                                              + u_prev_newton->dy[i] * v->dy[i])
                                + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] 
                                                                + u->dy[i] * v->dy[i])));
  return result;
}

// Residual vector.
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i]
                       + TAU * (lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] 
                                                              + u_prev_newton->dy[i] * v->dy[i])
		                - heat_src(e->x[i], e->y[i], 0) * v->val[i]));
  return result;
}

// Dummy bilinear and linar forms to get mass matrix
template<typename Real, typename Scalar>
Scalar dummy_bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (SIGMA * u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar dummy_linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (1);
  return result;
}


template<typename Real, typename Scalar>
Scalar jac_ss(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext, double t)
{
  Scalar result = 0;
  Func<Scalar>* Y1_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (BETA * u->val[i] * v->val[i]
                       + GAMMA * TAU * (dlam_du(Y1_prev_newton->val[i]) * u->val[i] * (Y1_prev_newton->dx[i] * v->dx[i]
                                                                                       + Y1_prev_newton->dy[i] * v->dy[i])
                                        + lam(Y1_prev_newton->val[i]) * (u->dx[i] * v->dx[i]
                                                                         + u->dy[i] * v->dy[i])));
  return result;
}

template<typename Real, typename Scalar>
Scalar res_ss(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext, double t)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i]
                                                      + u_prev_newton->dy[i] * v->dy[i])
                       - heat_src(e->x[i], e->y[i], t) * v->val[i]);
  return result;
}




//////////////////////////// SDIRK22 ////////////
// Y1 ///////////
template<typename Real, typename Scalar>
Scalar jac1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* Y1_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (BETA * u->val[i] * v->val[i]
                       + GAMMA * TAU * (dlam_du(Y1_prev_newton->val[i]) * u->val[i] * (Y1_prev_newton->dx[i] * v->dx[i] 
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
    result += wt[i] * (BETA * (Y1_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] 
                       + GAMMA * TAU * (lam(Y1_prev_newton->val[i]) * (Y1_prev_newton->dx[i] * v->dx[i] 
                                                                       + Y1_prev_newton->dy[i] * v->dy[i])
                                        - heat_src(e->x[i], e->y[i], 0) * v->val[i]));
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
                       + GAMMA * TAU * (dlam_du(Y2_prev_newton->val[i]) * u->val[i] * (Y2_prev_newton->dx[i] * v->dx[i]
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
                       + GAMMA * TAU * (lam(Y2_prev_newton->val[i]) * (Y2_prev_newton->dx[i] * v->dx[i]
                                                                 + Y2_prev_newton->dy[i] * v->dy[i])
                                        - heat_src(e->x[i], e->y[i], 0) * v->val[i]));
  return result;
}

