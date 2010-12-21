// Bilinear volume forms.
template<typename Real, typename Scalar>
Scalar bilinear_form_vol_SE(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_SE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
template<typename Real, typename Scalar>
Scalar bilinear_form_vol_NE(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_NE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
template<typename Real, typename Scalar>
Scalar bilinear_form_vol_NW(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_NW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
template<typename Real, typename Scalar>
Scalar bilinear_form_vol_SW(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_SW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }

// Linear volume forms.
template<typename Real, typename Scalar>
Scalar linear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return RHS * int_v<Real, Scalar>(n, wt, v); }

// Linear surface forms.
template<typename Real, typename Scalar>
Scalar linear_form_surf_VERTICAL_SE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return - A_SE * int_v<Real, Scalar>(n, wt, v); }
template<typename Real, typename Scalar>
Scalar linear_form_surf_VERTICAL_NE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return - A_NE * int_v<Real, Scalar>(n, wt, v); }
template<typename Real, typename Scalar>
Scalar linear_form_surf_VERTICAL_NW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return - A_NW * int_v<Real, Scalar>(n, wt, v); }
template<typename Real, typename Scalar>
Scalar linear_form_surf_VERTICAL_SW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return - A_SW * int_v<Real, Scalar>(n, wt, v); }
template<typename Real, typename Scalar>
Scalar linear_form_surf_TOP_NE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_NE * int_v<Real, Scalar>(n, wt, v); }
template<typename Real, typename Scalar>
Scalar linear_form_surf_TOP_NW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{ return A_NW * int_v<Real, Scalar>(n, wt, v); }