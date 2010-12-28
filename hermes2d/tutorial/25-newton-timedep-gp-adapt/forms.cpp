// Residual for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar F_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_iter = u_ext[0];
  Func<Scalar>* sln_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * (psi_iter->val[i] - sln_prev_time->val[i]) * v->val[i] / TAU
            - H*H/(2*M) * (psi_iter->dx[i] * v->dx[i] + psi_iter->dy[i] * v->dy[i])
            - G * psi_iter->val[i] *  psi_iter->val[i] * conj(psi_iter->val[i]) * v->val[i]
            - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * psi_iter->val[i] * v->val[i]);

  return result;
}

template<typename Real, typename Scalar>
Scalar J_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_iter = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * u->val[i] * v->val[i] / TAU
                     - H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                     - 2* G * u->val[i] *  psi_iter->val[i] * conj(psi_iter->val[i]) * v->val[i]
                     - G * psi_iter->val[i] * psi_iter->val[i] * u->val[i] * v->val[i]
                     - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar F_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_iter = u_ext[0];
  Func<Scalar>* sln_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * (psi_iter->val[i] - sln_prev_time->val[i]) * v->val[i] / TAU
            - 0.5*H*H/(2*M) * (psi_iter->dx[i] * v->dx[i] + psi_iter->dy[i] * v->dy[i])
            - 0.5*H*H/(2*M) * (sln_prev_time->dx[i] * v->dx[i] + sln_prev_time->dy[i] * v->dy[i])
            - 0.5*G * psi_iter->val[i] *  psi_iter->val[i] * conj(psi_iter->val[i]) * v->val[i]
            - 0.5*G * sln_prev_time->val[i] *  sln_prev_time->val[i] * conj(sln_prev_time->val[i]) * v->val[i]
            - 0.5*0.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * (psi_iter->val[i] + sln_prev_time->val[i]) * v->val[i]);

  return result;
}

template<typename Real, typename Scalar>
Scalar J_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_iter = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * u->val[i] * v->val[i] / TAU
                     - 0.5*H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                     - 0.5*2* G * u->val[i] *  psi_iter->val[i] * conj(psi_iter->val[i]) * v->val[i]
                     - 0.5*G * psi_iter->val[i] *  psi_iter->val[i] * u->val[i] * v->val[i]
                     - 0.5*.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
  return result;
}

