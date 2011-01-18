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
