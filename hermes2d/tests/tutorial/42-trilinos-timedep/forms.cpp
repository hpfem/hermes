template<typename Real, typename Scalar>
Scalar jacobian(int n, double *wt, Func<Real> *u_ext[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * RHO * vi->val[i] * vj->val[i] / TAU
                     + LAMBDA * (vi->dx[i] * vj->dx[i] + vi->dy[i] * vj->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar jacobian_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ALPHA * LAMBDA * vi->val[i] * vj->val[i];
  return result;
}

template<typename Real, typename Scalar>
Scalar residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (HEATCAP * RHO * (u_ext[0]->val[i] - ext->fn[0]->val[i]) * vj->val[i] / TAU
		       + LAMBDA * (u_ext[0]->dx[i] * vj->dx[i] + u_ext[0]->dy[i] * vj->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar residual_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vj, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ALPHA * LAMBDA * (u_ext[0]->val[i] - TEMP_EXT) * vj->val[i];
  return result;
}
