// Bilinear form for the implicit Euler time discretization, upper left block.
template<typename Real, typename Scalar>
Scalar biform_euler_0_0(int n, double *wt, Func<Scalar> *u_ext[], 
                        Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result = result + wt[i] * (u->val[i] / TAU) * v->val[i];
  return result;
}

// Bilinear form for the implicit Euler time discretization, upper right block.
template<typename Real, typename Scalar>
Scalar biform_euler_0_1(int n, double *wt, Func<Scalar> *u_ext[], 
                        Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result = result + wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

// Bilinear form for the implicit Euler time discretization, lower left block.
template<typename Real, typename Scalar>
Scalar biform_euler_1_0(int n, double *wt, Func<Scalar> *u_ext[], 
                        Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result = result - wt[i] * u->val[i] * v->val[i];
  return result;
}

// Bilinear form for the implicit Euler time discretization, lower left block.
template<typename Real, typename Scalar>
Scalar biform_euler_1_1(int n, double *wt, Func<Scalar> *u_ext[], 
                        Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result = result + wt[i] * u->val[i] / TAU * v->val[i];
  return result;
}

// Right-hand side for the implicit Euler time discretization, first component.
template<typename Real, typename Scalar>
Scalar liform_euler_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* phi_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result = result + wt[i] * (phi_prev_time->val[i] / TAU) * v->val[i];

  return result;
}

// Right-hand side for the implicit Euler time discretization, second component.
template<typename Real, typename Scalar>
Scalar liform_euler_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* psi_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result = result + wt[i] * (psi_prev_time->val[i] / TAU) * v->val[i];

  return result;
}


