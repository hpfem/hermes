template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{  
  if (PROB_PARAM == 0){
   ALPHA = 20;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 1){
   ALPHA = 1000;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 2){
   ALPHA = 1000;
   X_LOC = 1.5;
   Y_LOC = 0.25;
   R_ZERO = 0.92;
   }
else{
   ALPHA = 50;
   X_LOC = 0.5;
   Y_LOC = 0.5;
   R_ZERO = 0.25;
   }

  Real a = pow(x - X_LOC, 2);
  Real b = pow(y - Y_LOC, 2);
  Real c = sqrt(a + b);
  Real d = ((ALPHA*x - (ALPHA * X_LOC)) * (2*x - (2 * X_LOC)));
  Real e = ((ALPHA*y - (ALPHA * Y_LOC)) * (2*y - (2 * Y_LOC)));
  Real f = (pow(ALPHA*c - (ALPHA * R_ZERO), 2) + 1.0);
  Real g = (ALPHA * c - (ALPHA * R_ZERO));

  return ((ALPHA/(c * f)) - (d/(2 * pow(a + b, 1.5) * f)) - ((ALPHA * d * g)/((a + b) * pow(f, 2))) + 
         (ALPHA/(c * f)) - (e/(2 * pow(a + b, 1.5) * f)) - ((ALPHA * e * g)/((a + b) * pow(f, 2))));
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

template<typename Real, typename Scalar>
Scalar residual_estimator(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Scalar result = 0.;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( rhs(e->x[i], e->y[i]) - u->laplace[i] );
  
  return result * sqr(e->diam);  
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif  
}


template<typename Real, typename Scalar>
Scalar interface_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                           e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i))  );
  return result * e->diam / 24.;
}