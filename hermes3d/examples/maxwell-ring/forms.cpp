
double jv(double n, double x);

// Weak forms.
template<typename real, typename scalar>
scalar biform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  cplx ii = cplx(0.0, 1.0);
  return -1.0/(mu_r*mu_0) * hcurl_int_curl_u_curl_v<real, scalar>(n, wt, u, v, e)
    + ii * OMEGA * SIGMA * hcurl_int_u_v<real, scalar>(n, wt, u, v, e);
}

/// Integral \F \v
///
template<typename f_t, typename res_t>
res_t hcurl_int_J_v(int n, double *wt, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++) {
		result += wt[i] * (v->val0[i] * J[0] + v->val1[i] * J[1] + v->val2[i] * J[2]);
	}
	return result;
}

template<typename real, typename scalar>
scalar liform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  return hcurl_int_J_v<real, scalar>(n, wt, v, e);
}


