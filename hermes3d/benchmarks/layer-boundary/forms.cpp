template<typename real, typename scalar>
scalar biform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, 
              Geom<real> *e, ExtData<scalar> *ext)
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e) + K*K * int_u_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar liform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
              Geom<real> *e, ExtData<scalar> *ext)
{
  return K*K * int_u<real, scalar>(n, wt, u, e);
}
