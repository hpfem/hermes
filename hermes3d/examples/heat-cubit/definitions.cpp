template<typename real, typename scalar>
scalar bilinear_form1(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                      ExtData<scalar> *data)
{
	return 10 * int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar bilinear_form2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                      ExtData<scalar> *data)
{
	return 0.5 * int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename T>
T rhs(T x, T y, T z)
{
	return 4.0;
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<real> *e, ExtData<scalar> *data)
{
	return -int_F_v<real, scalar>(n, wt, rhs, u, e);
}
