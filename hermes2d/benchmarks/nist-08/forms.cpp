template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar val = 0;
  for (int i=0; i < n; i++) {
    Scalar x = e->x[i];
    Scalar y = e->y[i];
    Scalar r = sqrt(x*x + y*y);
    Scalar h = 1/(ALPHA + r);
    Scalar grad_u_grad_v = u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - pow(h, 4) * u->val[i] * v->val[i]);
  }

  return val;
}

template<typename Real>
Real rhs(Real x, Real y)
{
  return -sin(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4) + 2*cos(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(1.0/2.0))) + pow(x,2)*sin(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2))) + pow(y,2)*sin(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2))) - pow(x,2)*cos(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0))) - pow(y,2)*cos(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0))) - 2*pow(x,2)*cos(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2))) - 2*pow(y,2)*cos(1.0/(ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((ALPHA + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2)));
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

// integration order for the linear_form
Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}
