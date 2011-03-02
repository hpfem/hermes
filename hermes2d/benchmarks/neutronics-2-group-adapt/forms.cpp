//////////   Eq 1   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  int m = e->elem_marker - 1;
  return (D[m][0]) * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
         (Sr[m][0] - chi[m][0]*nSf[m][0] ) * int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_0_0(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) * int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform_surf_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 8*D[e->elem_marker-1][0]*int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_surf_0_0(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform_0_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  int m = e->elem_marker - 1;
  return -chi[m][1]*nSf[m][1] * int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_0_1(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar liform_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, Q1, v, e); 
}

// Integration order for the volumetric linear form.
// Set to the polynomial degree of the test function plus a sufficiently high addition to correctly integrate the non-polynomial, 
// but still relatively smooth source function for 1st group.
Ord liform_0_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(20 + v->val[0].get_order()); 
}

//////////   Eq 2   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  int m = e->elem_marker - 1;
  return (D[m][1]) * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
         (Sr[m][1]) * int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_1_1(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) * int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform_surf_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 8*D[e->elem_marker-1][1] * int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_surf_1_1(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform_1_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -Ss[e->elem_marker-1][1][0] * int_u_v<Real, Scalar>(n, wt, u, v);
}
template<>
Ord biform_1_0(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_u_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar liform_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	return int_F_v<Real, Scalar>(n, wt, Q2, v, e); 
}

// Integration order for the volumetric linear form.
// Set to the polynomial degree of the test function plus a sufficiently high addition to correctly integrate the non-polynomial, 
// highly oscillating source function for 2nd group.
Ord liform_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30 + v->val[0].get_order());
}
