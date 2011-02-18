// Bilinear and linear forms corresponding to simple linearization
// of the convective term.
template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_0_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) / RE + int_u_v<Real, Scalar>(n, wt, u, v) / TAU;
}

template<typename Real, typename Scalar>
Scalar simple_bilinear_form_unsym_0_0_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  Func<Scalar>* xvel_prev_time = ext->fn[0];
  Func<Scalar>* yvel_prev_time = ext->fn[1];
  result = int_w_nabla_u_v<Real, Scalar>(n, wt, xvel_prev_time, yvel_prev_time, u, v);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar simple_linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* vel_prev_time = ext->fn[0]; // this form is used with both velocity components
  return int_u_v<Real, Scalar>(n, wt, vel_prev_time, v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_0_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  result = - int_u_dvdx<Real, Scalar>(n, wt, p, v);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_unsym_1_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *p, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  result = - int_u_dvdy<Real, Scalar>(n, wt, p, v);
#endif
  return result;
}

// Bilinear and linear forms corresponding to the Newton's method.
template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  Func<Scalar>* xvel_prev_newton = u_ext[0];
  Func<Scalar>* yvel_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
                        * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_0_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  Func<Scalar>* xvel_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_1_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  Func<Scalar>* yvel_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_unsym_1_1(int n, double *wt, Func<Scalar> *u_ext[],  Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
#ifndef STOKES
  Func<Scalar>* xvel_prev_newton = u_ext[0];
  Func<Scalar>* yvel_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
                       * v->val[i] * yvel_prev_newton->dy[i]);
#endif
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xvel_prev_time = ext->fn[0];  
  Func<Scalar>* yvel_prev_time = ext->fn[1];
  Func<Scalar>* xvel_prev_newton = u_ext[0];  
  Func<Scalar>* yvel_prev_newton = u_ext[1];  
  Func<Scalar>* p_prev_newton = u_ext[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / TAU +
                       (xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / RE 
#ifndef STOKES
                       + (xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i] 
#endif
                       - (p_prev_newton->val[i] * v->dx[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xvel_prev_time = ext->fn[0];  
  Func<Scalar>* yvel_prev_time = ext->fn[1];
  Func<Scalar>* xvel_prev_newton = u_ext[0];  
  Func<Scalar>* yvel_prev_newton = u_ext[1];  
  Func<Scalar>* p_prev_newton = u_ext[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / TAU +
                       (yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / RE 
#ifndef STOKES
                       + (xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i] 
#endif
                       - (p_prev_newton->val[i] * v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_F_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* xvel_prev_newton = u_ext[0];  
  Func<Scalar>* yvel_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++)
    result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
  return result;
}
