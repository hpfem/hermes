#ifdef RHS2

double fnc(double x, double y, double z)
{
	return x*x + y*y + z*z;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *data) 
{
	scalar res = 0.0;
	// FIXME: 
//	for (int i = 0; i < n; i++)
//		res += wt[i] * (data->ext[0].fn[i] * v->fn[i]);
	return res;
}

#elif defined SYS

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
//	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return x;
//	return 3.0;
//	return (1 - x*x) * x*x * (1 - y*y) * y*y * (1 - z*z) * z*z;
}

// needed for calculation norms and used by visualizator
double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = -2 * x * (1 - y*y) * (1 - z*z);
//	dy = -2 * (1 - x*x) * y * (1 - z*z);
//	dz = -2 * (1 - x*x) * (1 - y*y) * z;
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = 2 * x * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	dy = 2 * x*x * (1 - x*x) * y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y*y * z*z * (1 - z*z);
//	dz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z*z;
	dx = 1;
	dy = 0;
	dz = 0;

	return u2(x, y, z);
}

//

BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

BCType bc_types_2(int marker)
{
	if (marker == 3) return BC_NATURAL;
	else return BC_ESSENTIAL;
}

scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z) {
	return u2(x, y, z);
}

template<typename real, typename scalar>
scalar bilinear_form_1_1(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                        ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar bilinear_form_1_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                        ExtData<scalar> *data)
{
//	return 0.0;
	return int_u_v<real, scalar>(n, wt, u, v, e);
}

//template<typename T>
//T f1(T x, T y, T z)
//{
//	T ddxx = -2 * (1 - y*y) * (1 - z*z);
//	T ddyy = -2 * (1 - x*x) * (1 - z*z);
//	T ddzz = -2 * (1 - x*x) * (1 - y*y);
//
//	return -(ddxx + ddyy + ddzz) + u2(x, y, z);
//}

template<typename real, typename scalar>
scalar linear_form_1(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<realt> *e, ExtData<scalar> *data)
{
//	return int_F_v<f_t, res_t>(n, wt, f1, u, e);
//	return -3.0 * int_u<f_t, res_t>(n, wt, u, e);
	scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * ((-6.0 + (e->x[i]) * u->fn[i]);
	return res;
//	int_u<f_t, scalar>(n, wt, u, e);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//template<typename f_t, typename res_t>
//res_t bilinear_form_2_1(int n, double *wt, Func<res_t> *u_ext[], Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e,
//                        ExtData<res_t> *data)
//{
//	return int_u_v<f_t, res_t>(n, wt, u, v, e);
//}

template<typename real, typename scalar>
scalar bilinear_form_2_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                        ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

//template<typename T>
//T f2(T x, T y, T z)
//{
//	T ddxx = 2 * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	T ddyy = 2 * x*x * (1 - x*x) * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * z*z * (1 - z*z);
//	T ddzz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z;
//
//	return -(ddxx + ddyy + ddzz) + u1(x, y, z);
//}

template<typename real, typename scalar>
scalar linear_form_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<real> *e, ExtData<scalar> *data)
{
//	return int_F_v<f_t, res_t>(n, wt, f2, u, e);
//	return -6.0 * int_u<f_t, res_t>(n, wt, u, e);
	return 0;
}

template<typename real, typename scalar>
scalar linear_form_2_surf(int np, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<real> *e, ExtData<scalar> *data)
{
	scalar result = 0;
	for (int i = 0; i < np; i++) {
		res_t dx = 1;
		res_t dy = 0;
		res_t dz = 0;

		result += wt[i] * (u->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return result;
}

#elif defined SYS3

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u3(T x, T y, T z)
{
	return x;
}

double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = -2 * x * (1 - y*y) * (1 - z*z);
	dy = -2 * (1 - x*x) * y * (1 - z*z);
	dz = -2 * (1 - x*x) * (1 - y*y) * z;

	return u2(x, y, z);
}

double exact_sln_fn_3(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 1;
	dy = 0;
	dz = 0;

	return u3(x, y, z);
}

//

BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

BCType bc_types_2(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z)
{
	return 0;
}

BCType bc_types_3(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_3(int ess_bdy_marker, double x, double y, double z)
{
	return x;
}

// 1. eqn ------------------------------------------------------------------------------------------

template<typename real, typename scalar>
scalar biform_1_1(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                 ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar biform_1_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                 ExtData<scalar> *data)
{
	return int_u_v<real, scalar>(n, wt, u, v, e);
}

template<typename T>
T rhs1(T x, T y, T z)
{
	return -6.0 + u2(x, y, z);
}

template<typename real, typename scalar>
scalar liform_1(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
	scalar res = 0.0;
	// FIXME: 
//	for (int i = 0; i < n; i++)
//		res += wt[i] * (rhs1(e->x[i], e->y[i], e->z[i]) * v->fn[i]);
	return res;
}

// 2. eqn ------------------------------------------------------------------------------------------

template<typename real, typename scalar>
scalar biform_2_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                 ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar biform_2_3(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                 ExtData<scalar> *data)
{
	return int_u_v<real, scalar>(n, wt, u, v, e);
}

template<typename T>
T rhs2(T x, T y, T z)
{
	T laplace = 2 * (1 - y*y) * (1 - z*z) + 2 * (1 - x*x) * (1 - z*z) + 2 * (1 - x*x) * (1 - y*y);
	return laplace + u3(x, y, z);
}

template<typename real, typename scalar>
scalar liform_2(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
	return int_F_v<real, scalar>(n, wt, rhs2, v, e);
}

// 3. eqn ------------------------------------------------------------------------------------------

template<typename real, typename scalar>
scalar biform_3_3(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e,
                        ExtData<scalar> *data)
{
	return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}


#endif
