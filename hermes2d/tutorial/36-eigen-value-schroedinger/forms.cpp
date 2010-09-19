template<typename Real, typename Scalar>
Scalar bilinear_form_H(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}



template<typename Real, typename Scalar>
Scalar bilinear_form_V(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
		       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result +=pot(x, y)*u->val[i]*v->val[i]* wt[i];
  }
  return result;
}

// Integration order for the bilinear form
Ord bilinear_form_V_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}




template<typename Real, typename Scalar>
Scalar bilinear_form_U(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}
