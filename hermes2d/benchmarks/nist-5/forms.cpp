template<typename Real, typename Scalar>
Scalar biform1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (P_1 * u->dx[i] * v->dx[i] + Q_1 * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar biform2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (P_2 * u->dx[i] * v->dx[i] + Q_2 * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar biform3(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (P_3 * u->dx[i] * v->dx[i] + Q_3 * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar biform4(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (P_4 * u->dx[i] * v->dx[i] + Q_4 * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar biform5(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (P_5 * u->dx[i] * v->dx[i] + Q_5 * u->dy[i] * v->dy[i]);
  return result;
}

// 
//
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return F_1 * int_v<Real, Scalar>(n, wt, v);
} 
template<typename Real, typename Scalar>
Scalar linear_form_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return F_2 * int_v<Real, Scalar>(n, wt, v);
} 
template<typename Real, typename Scalar>
Scalar linear_form_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return F_3 * int_v<Real, Scalar>(n, wt, v);
} 
template<typename Real, typename Scalar>
Scalar linear_form_4(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return F_4 * int_v<Real, Scalar>(n, wt, v);
} 
template<typename Real, typename Scalar>
Scalar linear_form_5(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return F_5 * int_v<Real, Scalar>(n, wt, v);
} 

