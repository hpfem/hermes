// Residuum for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar F_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  Func<Scalar>* psi_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * (psi_prev_newton->val[i] - psi_prev_time->val[i]) * v->val[i] / time_step
            - H*H/(2*M) * (psi_prev_newton->dx[i] * v->dx[i] + psi_prev_newton->dy[i] * v->dy[i])
            - G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
            - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * psi_prev_newton->val[i] * v->val[i]);

  return result;
}

// Jacobian for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar J_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * u->val[i] * v->val[i] / time_step
                     - H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                     - 2* G * u->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                     - G * psi_prev_newton->val[i] * psi_prev_newton->val[i] * u->val[i] * v->val[i]
                     - .5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
  return result;
}

// Residuum for the Crank-Nicolson method
template<typename Real, typename Scalar>
Scalar F_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  Func<Scalar>* psi_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * (psi_prev_newton->val[i] - psi_prev_time->val[i]) * v->val[i] / time_step
            - 0.5*H*H/(2*M) * (psi_prev_newton->dx[i] * v->dx[i] + psi_prev_newton->dy[i] * v->dy[i])
            - 0.5*H*H/(2*M) * (psi_prev_time->dx[i] * v->dx[i] + psi_prev_time->dy[i] * v->dy[i])
            - 0.5*G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
            - 0.5*G * psi_prev_time->val[i] *  psi_prev_time->val[i] * conj(psi_prev_time->val[i]) * v->val[i]
            - 0.5*0.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * (psi_prev_newton->val[i] + psi_prev_time->val[i]) * v->val[i]);

  return result;
}

// Jacobian for the Crank-Nicolson method
template<typename Real, typename Scalar>
Scalar J_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (ii * H * u->val[i] * v->val[i] / time_step
                     - 0.5*H*H/(2*M) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                     - 0.5*2* G * u->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                     - 0.5*G * psi_prev_newton->val[i] *  psi_prev_newton->val[i] * u->val[i] * v->val[i]
                     - 0.5*.5*M*OMEGA*OMEGA * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
  return result;
}

