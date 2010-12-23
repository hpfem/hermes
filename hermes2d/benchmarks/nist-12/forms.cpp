template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  //For more elegant showing please execute file "generate_rhs.py" 

  Real a_P = (-ALPHA_P * pow((x - X_P), 2) - ALPHA_P * pow((y - Y_P), 2));

  Real a_W = pow(x - X_W, 2);
  Real b_W = pow(y - Y_W, 2);
  Real c_W = sqrt(a_W + b_W);
  Real d_W = ((ALPHA_W * x - (ALPHA_W * X_W)) * (2 * x - (2 * X_W)));
  Real e_W = ((ALPHA_W * y - (ALPHA_W * Y_W)) * (2 * y - (2 * Y_W)));
  Real f_W = (pow(ALPHA_W * c_W - (ALPHA_W * R_0), 2) + 1.0);
  Real g_W = (ALPHA_W * c_W - (ALPHA_W * R_0));

  return 4 * exp(a_P) * ALPHA_P * (ALPHA_P * (x - X_P) * (x - X_P) + ALPHA_P * (y - Y_P) * (y - Y_P) - 1)
         + ((ALPHA_W/(c_W * f_W)) - (d_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((ALPHA_W * d_W * g_W)/((a_W + b_W) * pow(f_W, 2))) 
         + (ALPHA_W/(c_W * f_W)) - (e_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((ALPHA_W * e_W * g_W)/((a_W + b_W) * pow(f_W, 2))))
         + (1.0 / EPSILON) * (1.0 / EPSILON) * exp(-(1 + y) / EPSILON);  
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
