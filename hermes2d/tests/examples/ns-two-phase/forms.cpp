/* THIS IS THE BASIS FOR A NEW WORKING CODE, REPLACING THE 
   FIRST UNCOMMENTED FUNCTION BELOW
double int_sign_u_v(double a, double b, Func<double>* l, Func<double>* fu, Func<double>* fv)
{
   Quad2D* quad = fu->get_quad_2d();
   int o = 5;
   limit_order(o);
   l->set_quad_order(o, H2D_FN_VAL);
   fu->set_quad_order(o);
   fv->set_quad_order(o);
   double* uval = fu->get_fn_values();
   double* vval = fv->get_fn_values();
   double* lval = l->get_fn_values();

  double result = 0;
  h1_integrate_expression(((lval[i] < 0.0) ? b : a) * uval[i] * vval[i]);
  return result;
}

Ord int_sign_u_v(double a, double b, Func<Ord>* l, Func<Ord>* fu, Func<Ord>* fv)
{
   Quad2D* quad = fu->get_quad_2d();
   int o = 5;
   limit_order(o);
   l->set_quad_order(o, H2D_FN_VAL);
   fu->set_quad_order(o);
   fv->set_quad_order(o);
   Ord* uval = fu->get_fn_values();
   Ord* vval = fv->get_fn_values();
   Ord* lval = l->get_fn_values();

   Ord result =0.0;
   h1_integrate_expression(uval[i] * vval[i]);  // For the integration it does not matter whether 
                                                // there is 'a' or 'b'.
  return result;
}
*/

template<typename Real, typename Scalar>
Scalar int_sign_u_v(int n, double *wt, double a, double b, 
                    Func<Real>* l, Func<Real>* fu, Func<Real>* fv)
{
  //h1_integrate_expression(((lval[i] < 0.0) ? b : a) * uval[i] * vval[i]);
  Scalar result = 0.0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (((l->val[i] < 0.0) ? b : a) * fu->val[i] * fv->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar int_sign_grad_u_grad_v(int n, double *wt, double nu1, double nu2, 
                              Func<Real>* l, Func<Real>* fu, Func<Real>* fv)
{
  //h1_integrate_dd_expression(((lval[i] < 0.0) ? nu2 : nu1) * (t_dudx * t_dvdx + t_dudy * t_dvdy));
  Scalar result = 0;
  for (int i = 0; i < n; i++)
   {
     result += wt[i] * (((l->val[i] < 0.0) ? nu2 : nu1) * (fu->dx[i] * fv->dx[i] + fu->dy[i] * fv->dy[i]));
   }
  return result;
}

template<typename Real, typename Scalar>
Scalar int_sign_w_nabla_u_v(int n, double *wt, double ro1, double ro2, 
                            Func<Real>* l, Func<Real>* w1, Func<Real>* w2, Func<Real>* fu, Func<Real>* fv)
{
  //h1_integrate_dd_expression(((lval[i] < 0.0) ? ro2 : ro1) * (w1val[i] * t_dudx + w2val[i] * t_dudy) * vval[i]);
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (((l->val[i] < 0.0) ? ro2 : ro1) * (w1->val[i] * fu->dx[i] + w2->val[i] * fu->dy[i]) * fv->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar int_sign_v(int n, double *wt, double a, double b, Func<Real>* l, Func<Real>* fu)
{
  //h1_integrate_expression(((lval[i] < 0.0) ? b : a) * uval[i]);)
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (((l->val[i] < 0.0) ? b : a) * fu->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real>* fu, Func<Real>* fv, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* xprev = ext->fn[1];
  Func<Scalar>* yprev = ext->fn[2];

  return int_sign_grad_u_grad_v<Real, Scalar>(n, wt, Nu1, Nu2, lprev, fu, fv) + 
         int_sign_u_v<Real, Scalar>(n, wt, Ro1, Ro2, lprev, fu, fv)/tau + 
         int_sign_w_nabla_u_v<Real, Scalar>(n, wt, Ro1, Ro2, lprev, xprev, yprev, fu, fv); 
}

//scalar bilinear_form_0_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return int_dudx_v(fu, fv, ru ,rv ); }
template<typename Real, typename Scalar>
Scalar bilinear_form_0_2(int n, double *wt, Func<Scalar> *u_ext[], 
                         Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return int_dudx_v<Real, Scalar>(n, wt, fu, fv); 
}
 
//scalar bilinear_form_1_2(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return int_dudy_v(fu, fv, ru ,rv ); } 
template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Scalar> *u_ext[], 
                         Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  return int_dudy_v<Real, Scalar>(n, wt, fu, fv); 
}

//scalar bilinear_form_2_0(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv) 
//  { return int_dudx_v(fu, fv, ru, rv);}
template<typename Real, typename Scalar>
Scalar bilinear_form_2_0(int n, double *wt, Func<Scalar> *u_ext[], 
                         Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext) 
{ 
  return int_dudx_v<Real, Scalar>(n, wt, fu, fv);
}

//scalar bilinear_form_2_1(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)  
//  { return int_dudy_v(fu, fv, ru, rv); }
template<typename Real, typename Scalar>
Scalar bilinear_form_2_1(int n, double *wt, Func<Scalar> *u_ext[], 
                         Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)  
{ 
  return int_dudy_v<Real, Scalar>(n, wt, fu, fv); 
}

//scalar bilinear_form_3_3(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
//  { return  int_u_v(fu, fv, ru, rv)/tau + int_w_nabla_u_v(&xprev, &yprev, fu, fv, ru, rv); }
template<typename Real, typename Scalar>
Scalar bilinear_form_3_3(int n, double *wt, Func<Scalar> *u_ext[], 
                         Func<Real>* fu, Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* xprev = ext->fn[0];
  Func<Scalar>* yprev = ext->fn[1];

  return  int_u_v<Real, Scalar>(n, wt, fu, fv)/tau 
          + int_w_nabla_u_v<Real, Scalar>(n, wt, xprev, yprev, fu, fv); 
}

//scalar linear_form_0(RealFunction* fv, RefMap* rv)
//  { return int_sign_u_v(Ro1, Ro2, &lprev, &xprev, fv, rv, rv)/tau; }
template<typename Real, typename Scalar>
Scalar linear_form_0(int n, double *wt, Func<Scalar> *u_ext[], 
                     Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* xprev = ext->fn[1];

  return int_sign_u_v<Real, Scalar>(n, wt, Ro1, Ro2, lprev, xprev, fv)/tau; 
}

//scalar linear_form_1(RealFunction* fv, RefMap* rv)
//  { return int_sign_u_v(Ro1, Ro2, &lprev, &yprev, fv, rv, rv)/tau - int_sign_v(10*Ro1, 10*Ro2, &lprev, fv, rv); }
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Scalar> *u_ext[], 
                     Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{ 
  Func<Scalar>* lprev = ext->fn[0];
  Func<Scalar>* yprev = ext->fn[1];

  return int_sign_u_v<Real, Scalar>(n, wt, Ro1, Ro2, lprev, yprev, fv)/tau 
         - int_sign_v<Real, Scalar>(n, wt, 10*Ro1, 10*Ro2, lprev, fv); 
}

//Scalar linear_form_3(RealFunction* fv, RefMap* rv)
//  { return int_u_v(&lprev, fv, rv, rv)/tau; }
template<typename Real, typename Scalar>
Scalar linear_form_3(int n, double *wt, Func<Scalar> *u_ext[], 
                     Func<Real>* fv, Geom<Real> *e, ExtData<Scalar> *ext)
{  
  Func<Scalar>* lprev = ext->fn[0];

  return int_u_v<Real, Scalar>(n, wt, lprev, fv)/tau; 
}
