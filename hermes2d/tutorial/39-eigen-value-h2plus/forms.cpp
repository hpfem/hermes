template<typename Real, typename Scalar>
Scalar bilinear_form_H(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)

{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += 4 * M_PI * x * wfun(x,y) * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) * wt[i];
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_V(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
		       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += pot(x, y) * 4 * M_PI * x * wfun(x,y) * u->val[i] * v->val[i]* wt[i];
  }
  return result;
}

// Integration order for the bilinear form.
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  // Returning the sum of the degrees of the basis and test function plus six.
  return u->val[0] * v->val[0] * e->x[0] * e->x[0] * e->x[0] * e->x[0] * e->x[0] * e->x[0]; 
}

template<typename Real, typename Scalar>
Scalar bilinear_form_U(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
		       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += 4 * M_PI * x * wfun(x,y) * u->val[i] * v->val[i] * wt[i];
  }
  return result;
}



