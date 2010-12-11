// Weak forms
template<typename Real, typename Scalar>
Scalar bilinear_form_surf_left(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double x = e->x[i];
    double y = e->y[i];
    double P = 0.0;
    double Q = 0.0;

    if (x == 0.0)
    {
      if ((y >= 0.0 && y <= 0.8)||(y >= 23.2 && y <= 24.0))
      {
        P = P_1; 
        Q = Q_1;
      }
      if ((y >= 1.6 && y <= 3.6)||(y >= 18.8 && y <= 21.2))
      {
        P = P_2; 
        Q = Q_2;
      }
      if (y >= 3.6 && y <= 18.8)
      {
        P = P_3; 
        Q = Q_3;
      }
      if ((y >= 0.8 && y <= 1.6)||(y >= 21.2 && y <= 23.2))
      {
        P = P_5; 
        Q = Q_5;
      }
    }
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_LEFT * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_right(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                                Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_RIGHT * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                              Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_TOP * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                                 Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_BOTTOM * u->val[i] * v->val[i]);
  }
  return result;
}

Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
//return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_left(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_LEFT * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_right(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_RIGHT * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_TOP * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_BOTTOM * int_v<Real, Scalar>(n, wt, v);
}

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

