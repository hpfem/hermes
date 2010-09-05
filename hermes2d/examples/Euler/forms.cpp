#include "hermes2d.h"
#include "params.h"

// Forward declaration.
double f_x(int i, double w0, double w1, double w3, double w4);

// Ondrej's code starts here.
			double matrix_R(int i, int j, double w0, double w1, double w3, double w4)
			{
					double rho = w0;
					double u = w1/w0;
					double w = w3/w0;
					double E = w4;
					double v2 = u*u+w*w;
					double p = (kappa-1)*(E - rho*v2/2);
					double c = sqrt(kappa*p/rho);
					if (i == 0 && j == 0)
							return 1;
					else if (i == 0 && j == 1)
							return 1;
					else if (i == 0 && j == 2)
							return 1;
					else if (i == 0 && j == 3)
							return 1;

					else if (i == 1 && j == 0)
							return u-c;
					else if (i == 1 && j == 1)
							return u;
					else if (i == 1 && j == 2)
							return u;
					else if (i == 1 && j == 3)
							return u+c;

					else if (i == 2 && j == 0)
							return w;
					else if (i == 2 && j == 1)
							return w;
					else if (i == 2 && j == 2)
							return w-c;
					else if (i == 2 && j == 3)
							return w;

					else if (i == 3 && j == 0)
							return v2/2 + c*c/(kappa-1) - u*c;
					else if (i == 3 && j == 1)
							return v2/2;
					else if (i == 3 && j == 2)
							return v2/2 - w*c;
					else if (i == 3 && j == 3)
							return v2/2 + c*c/(kappa-1) + u*c;

					printf("i=%d, j=%d;\n", i, j);
					error("Invalid index.");
			}

			double matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4)
			{
					double rho = w0;
					double u = w1/w0;
					double w = w3/w0;
					double E = w4;
					double v2 = u*u+w*w;
					double p = (kappa-1)*(E - rho*v2/2);
					double c = sqrt(kappa*p/rho);
					double result;
					if (i == 0 && j == 0)
							result = ((kappa-1)*v2/2 + u*c)/2;
					else if (i == 0 && j == 1)
							result = -(c+u*(kappa-1))/2;
					else if (i == 0 && j == 2)
							result = -w*(kappa-1)/2;
					else if (i == 0 && j == 3)
							result = (kappa-1)/2;

					else if (i == 1 && j == 0)
							result = c*c-c*w-(kappa-1)*v2/2;
					else if (i == 1 && j == 1)
							result = u*(kappa-1);
					else if (i == 1 && j == 2)
							result = c+w*(kappa-1);
					else if (i == 1 && j == 3)
							result = 1-kappa;

					else if (i == 2 && j == 0)
							result = w*c;
					else if (i == 2 && j == 1)
							result = 0;
					else if (i == 2 && j == 2)
							result = -c;
					else if (i == 2 && j == 3)
							result = 0;

					else if (i == 3 && j == 0)
							result = ((kappa-1)*v2/2 - u*c)/2;
					else if (i == 3 && j == 1)
							result = (c-u*(kappa-1))/2;
					else if (i == 3 && j == 2)
							result = -w*(kappa-1)/2;
					else if (i == 3 && j == 3)
							result = (kappa-1)/2;
					else {
							printf("i=%d, j=%d;\n", i, j);
							error("Invalid index.");
					}
					return result/(c*c);
			}

			double matrix_D_minus(int i, int j, double w0, double w1, double w3, double w4)
			{
					double rho = w0;
					double u = w1/w0;
					double w = w3/w0;
					double E = w4;
					double v2 = u*u+w*w;
					double p = (kappa-1)*(E - rho*v2/2);
					double c = sqrt(kappa*p/rho);
					double u_diag = 0;
					if (u < 0)
							u_diag = u;
					if (i == 0 && j == 0)
							return u-c;
					else if (i == 0 && j == 1)
							return 0;
					else if (i == 0 && j == 2)
							return 0;
					else if (i == 0 && j == 3)
							return 0;

					else if (i == 1 && j == 0)
							return 0;
					else if (i == 1 && j == 1)
							return u_diag;
					else if (i == 1 && j == 2)
							return 0;
					else if (i == 1 && j == 3)
							return 0;

					else if (i == 2 && j == 0)
							return 0;
					else if (i == 2 && j == 1)
							return 0;
					else if (i == 2 && j == 2)
							return u_diag;
					else if (i == 2 && j == 3)
							return 0;

					else if (i == 3 && j == 0)
							return 0;
					else if (i == 3 && j == 1)
							return 0;
					else if (i == 3 && j == 2)
							return 0;
					else if (i == 3 && j == 3)
							return 0;

					printf("i=%d, j=%d;\n", i, j);
					error("Invalid index.");
			}

			// multiplies two matrices
			void dot(double result[4][4], double A[4][4], double B[4][4])
			{
					for (int i=0; i < 4; i++)
							for (int j=0; j < 4; j++) {
									double sum=0;
									for (int k=0; k < 4; k++)
											sum += A[i][k] * B[k][j];
									result[i][j] = sum;
							}
			}

			// multiplies a matrix and a vector
			void dot_vector(double result[4], double A[4][4], double B[4])
			{
					for (int i=0; i < 4; i++) {
							double sum=0;
							for (int k=0; k < 4; k++)
									sum += A[i][k] * B[k];
							result[i] = sum;
					}
			}


			void T_rot(double result[4][4], double beta)
			{
					for (int i=0; i < 4; i++)
							for (int j=0; j < 4; j++)
									result[i][j] = 0;
					result[0][0] = 1;
					result[1][1] = cos(beta);
					result[1][2] = sin(beta);
					result[2][1] = -sin(beta);
					result[2][2] = cos(beta);
					result[3][3] = 1;
			}

			void A_minus(double result[4][4], double w0, double w1, double w3, double w4)
			{
					double _R[4][4];
					double _D_minus[4][4];
					double _R_inv[4][4];
					double _tmp[4][4];
					for (int i=0; i < 4; i++)
							for (int j=0; j < 4; j++)
									_R[i][j] = matrix_R(i, j, w0, w1, w3, w4);
					for (int i=0; i < 4; i++)
							for (int j=0; j < 4; j++)
									_D_minus[i][j] = matrix_D_minus(i, j, w0, w1, w3, w4);
					for (int i=0; i < 4; i++)
							for (int j=0; j < 4; j++)
									_R_inv[i][j] = matrix_R_inv(i, j, w0, w1, w3, w4);
					dot(_tmp, _D_minus, _R_inv);
					dot(result, _R, _tmp);
			}

			void riemann_solver(double result[4], double w_l[4], double w_r[4])
			{
					double _tmp1[4][4];
					double _tmp2[4][4];
					double _tmp3[4];
					double _tmp4[4];
					A_minus(_tmp1, w_r[0], w_r[1], w_r[2], w_r[3]);
					A_minus(_tmp2, w_l[0], w_l[1], w_l[2], w_l[3]);
					dot_vector(_tmp3, _tmp1, w_r);
					dot_vector(_tmp4, _tmp2, w_l);
					for (int i=0; i < 4; i++) 
					{
							double _1 = f_x(i, w_l[0], w_l[1], w_l[2], w_l[3]);
							double _2 = _tmp3[i];
							double _3 = _tmp4[i];
							result[i] = _1 + _2 - _3;
					}
			}

			// calculates the inverted flux, for testing purposes
			// it should return the same thing as riemann_solver(), only with minus sign
			void riemann_solver_invert(double result[4], double w_l[4], double w_r[4])
			{
					double m[4][4];
					double _w_l[4];
					double _w_r[4];
					double _tmp[4];
					T_rot(m, M_PI);
					dot_vector(_w_l, m, w_l);
					dot_vector(_w_r, m, w_r);
					riemann_solver(_tmp, _w_r, _w_l);
					T_rot(m, -M_PI);
					dot_vector(result, m, _tmp);
			}

			// Calculates the numerical flux in the normal (nx, ny) by rotating into the
			// local system, solving the Riemann problem and rotating back. It returns the
			// state.
			double numerical_flux(int i, double w_l[4], double w_r[4], double nx, double ny)
			{
				double result[4];
				double alpha = atan2(ny, nx);
				double mat_rot[4][4];
				double mat_rot_inv[4][4];
				double w_l_local[4];
				double w_r_local[4];
				double flux_local[4];
				T_rot(mat_rot, alpha);
				T_rot(mat_rot_inv, -alpha);
				dot_vector(w_l_local, mat_rot, w_l);
				dot_vector(w_r_local, mat_rot, w_r);
				riemann_solver(flux_local, w_l_local, w_r_local);
				dot_vector(result, mat_rot_inv, flux_local);
				return result[i];
			}

// Ondrej's code ends here.

///////////////////////////////////////////////////
////////Time dependent boundary conditions/////////
///////////////////////////////////////////////////
double bc_density(double y)
{
	return 1.0;
}

double bc_density_vel1(double y)
{
	return 50.0 * (1-y) * y;
}

double bc_density_vel2(double y)
{
	return 0;
}

double bc_energy(double y)
{
	return 1E5 + 0.5 * bc_density(y) * bc_density_vel1(y) * bc_density_vel1(y);
}

double bc_pressure(double y)
{
	return R/c_v * (bc_energy(y) - (bc_density_vel1(y) * bc_density_vel1(y) + bc_density_vel2(y) * bc_density_vel2(y))/(2 * bc_density(y)));
}

double ic_density(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density(y);
}
double ic_density_vel1(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density_vel1(y);
}
double ic_density_vel2(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density_vel2(y);
}
double ic_energy(double x, double y, scalar& dx, scalar& dy)
{
	return bc_energy(y);
}


template<typename Scalar>
Scalar calc_pressure(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return (R/c_v) * (e - (d_v1*d_v1 + d_v2*d_v2) / (2*d));
}

template<typename Scalar>
Scalar calc_a(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return std::sqrt((1+R/c_v) * calc_pressure(d,d_v1,d_v2,e) / d);
}

template<typename Scalar>
Scalar calc_energy(Scalar d, Scalar d_v1, Scalar d_v2, Scalar p)
{
	return (c_v/R) * p + (d_v1*d_v1+d_v2*d_v2) / 2*d;
}

///////////////////////////////////////////////////
////////First flux jacobian/////////////////////
///////////////////////////////////////////////////
template<typename Scalar>
Scalar A_1_0_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_0_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 1;
}

template<typename Scalar>
Scalar A_1_0_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_0_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_1_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - ((d_v1 * d_v1) / (d * d)) + 0.5 * (R / c_v) * ((d_v1 * d_v1 + d_v2 * d_v2) / 	(d * d));
}

template<typename Scalar>
Scalar A_1_1_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 2 * (d_v1 / d) - (R / c_v) * (d_v1 / d);
}

template<typename Scalar>
Scalar A_1_1_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (R / c_v) * (d_v2 / d);;
}

template<typename Scalar>
Scalar A_1_1_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return (R / c_v);
}

template<typename Scalar>
Scalar A_1_2_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - d_v1 * d_v2 / (d * d);
}

template<typename Scalar>
Scalar A_1_2_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v2 / d;
}

template<typename Scalar>
Scalar A_1_2_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v1 / d;
}

template<typename Scalar>
Scalar A_1_2_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_3_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (d_v1 * e) / (d * d) - (d_v1 / (d * d)) * (R / c_v) * (e - ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d))) + (d_v1 / d) * (R / c_v) * ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d * d));
}

template<typename Scalar>
Scalar A_1_3_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return (e / d) + (1 / d) * (R / c_v) * ( e - ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d * d))) - (R / c_v) * ((d_v1 * d_v1) / (d * d));
}

template<typename Scalar>
Scalar A_1_3_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (R / c_v) * (d_v1 * d_v2) / (d * d);
}

template<typename Scalar>
Scalar A_1_3_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v1 / d + (R / c_v) * (d_v1 / d);
}


///////////////////////////////////////////////////
////////First flux////////////////////////////////
///////////////////////////////////////////////////

template<typename Scalar>
Scalar f_1_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return A_1_0_0<Scalar>(d, d_v1, d_v2, e) * d 
		+ A_1_0_1<Scalar>(d, d_v1, d_v2, e) * d_v1 
		+ A_1_0_2<Scalar>(d, d_v1, d_v2, e) * d_v2
		+ A_1_0_3<Scalar>(d, d_v1, d_v2, e) * e;
}

template<typename Scalar>
Scalar f_1_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return A_1_1_0<Scalar>(d, d_v1, d_v2, e) * d 
		+ A_1_1_1<Scalar>(d, d_v1, d_v2, e) * d_v1 
		+ A_1_1_2<Scalar>(d, d_v1, d_v2, e) * d_v2
		+ A_1_1_3<Scalar>(d, d_v1, d_v2, e) * e;
}

template<typename Scalar>
Scalar f_1_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return A_1_2_0<Scalar>(d, d_v1, d_v2, e) * d 
		+ A_1_2_1<Scalar>(d, d_v1, d_v2, e) * d_v1 
		+ A_1_2_2<Scalar>(d, d_v1, d_v2, e) * d_v2
		+ A_1_2_3<Scalar>(d, d_v1, d_v2, e) * e;
}

template<typename Scalar>
Scalar f_1_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return A_1_3_0<Scalar>(d, d_v1, d_v2, e) * d 
		+ A_1_3_1<Scalar>(d, d_v1, d_v2, e) * d_v1 
		+ A_1_3_2<Scalar>(d, d_v1, d_v2, e) * d_v2
		+ A_1_3_3<Scalar>(d, d_v1, d_v2, e) * e;
}

double f_x(int i, double w0, double w1, double w3, double w4)
{
	if(i == 0)
		return f_1_0<double>(w0,w1,w3,w4);
	if(i == 1)
		return f_1_1<double>(w0,w1,w3,w4);
	if(i == 2)
		return f_1_2<double>(w0,w1,w3,w4);
	if(i == 3)
		return f_1_3<double>(w0,w1,w3,w4);
}

///////////////////////////////////////////////////
////////Second flux jacobian/////////////////////
///////////////////////////////////////////////////


template<typename Scalar>
Scalar A_2_0_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_0_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_0_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 1;
}

template<typename Scalar>
Scalar A_2_0_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_1_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - d_v1 * d_v2 / (d * d);
}

template<typename Scalar>
Scalar A_2_1_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v2 / d;
}

template<typename Scalar>
Scalar A_2_1_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v1 / d;
}

template<typename Scalar>
Scalar A_2_1_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_2_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - ((d_v2 * d_v2) / (d * d)) + 0.5 * (R / c_v) * ((d_v1 * d_v1 + d_v2 * d_v2) / 	(d * d));

}

template<typename Scalar>
Scalar A_2_2_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (R / c_v) * (d_v1 / d);
	
}

template<typename Scalar>
Scalar A_2_2_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return 2 * (d_v2 / d) - (R / c_v) * (d_v2 / d);
}

template<typename Scalar>
Scalar A_2_2_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return R / c_v;
}

template<typename Scalar>
Scalar A_2_3_0(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (d_v2 * e) / (d * d) - (d_v2 / (d * d)) * (R / c_v) * (e - ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d))) + (d_v2 / d) * (R / c_v) * ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d * d));
}

template<typename Scalar>
Scalar A_2_3_1(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return - (R / c_v) * (d_v1 * d_v2) / (d * d);
}

template<typename Scalar>
Scalar A_2_3_2(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return (e / d) + (1 / d) * (R / c_v) * ( e - ((d_v1 * d_v1 + d_v2 * d_v2) / (2 * d * d))) - (R / c_v) * ((d_v2 * d_v2) / (d * d));
}

template<typename Scalar>
Scalar A_2_3_3(Scalar d, Scalar d_v1, Scalar d_v2, Scalar e)
{
	return d_v2 / d + (R / c_v) * (d_v2 / d);
}


///////////////////////////////////////////////////
////////Volumetric bilinear forms//////////////////
///////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->dx[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_2(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->dy[i]);
  return result;
}

/*
template<typename Real, typename Scalar>
Scalar bilinear_form_0_3(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * ( 
  return result;
}
*/

template<typename Real, typename Scalar>
Scalar bilinear_form_1_0(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i] / TAU;

	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_3(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_0(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_1(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_2(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i] / TAU;

	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_3(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_3_0(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i] + A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_3_1(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i] + A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_3_2(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_3_3(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i] / TAU;

	for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i] + A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i]);

  return result;
}


///////////////////////////////////////////////////
////////Surface bilinear forms//////////////////
///////////////////////////////////////////////////

Ord linear_form_order(int n, double *wt, Func<Ord> *ue[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
	return Ord(20);
}

double linear_form_interface_0(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w02, w11, w12, w21, w22, w31, w32;
	double w_l[4], w_r[4];
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->get_val_central(i);
		w02 = ext->fn[0]->get_val_neighbor(i);
		
		w11 = ext->fn[1]->get_val_central(i);
		w12 = ext->fn[1]->get_val_neighbor(i);

		w21 = ext->fn[2]->get_val_central(i);
		w22 = ext->fn[2]->get_val_neighbor(i);

		w31 = ext->fn[3]->get_val_central(i);
		w32 = ext->fn[3]->get_val_neighbor(i);

		w_l[0] = w01;
		w_l[1] = w11;
		w_l[2] = w21;
		w_l[3] = w31;

		w_r[0] = w02;
		w_r[1] = w12;
		w_r[2] = w22;
		w_r[3] = w32;

		result += wt[i] * v->val[i] * numerical_flux(0,w_l,w_r,e->nx[i], e->ny[i]);
  }
  return result;
}

double linear_form_interface_1(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w02, w11, w12, w21, w22, w31, w32;
	double w_l[4], w_r[4];
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->get_val_central(i);
		w02 = ext->fn[0]->get_val_neighbor(i);
		
		w11 = ext->fn[1]->get_val_central(i);
		w12 = ext->fn[1]->get_val_neighbor(i);

		w21 = ext->fn[2]->get_val_central(i);
		w22 = ext->fn[2]->get_val_neighbor(i);

		w31 = ext->fn[3]->get_val_central(i);
		w32 = ext->fn[3]->get_val_neighbor(i);

		w_l[0] = w01;
		w_l[1] = w11;
		w_l[2] = w21;
		w_l[3] = w31;

		w_r[0] = w02;
		w_r[1] = w12;
		w_r[2] = w22;
		w_r[3] = w32;

		result += wt[i] * v->val[i] * numerical_flux(1,w_l,w_r,e->nx[i], e->ny[i]);
  }
  return result;
}

double linear_form_interface_2(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w02, w11, w12, w21, w22, w31, w32;
	double w_l[4], w_r[4];
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->get_val_central(i);
		w02 = ext->fn[0]->get_val_neighbor(i);
		
		w11 = ext->fn[1]->get_val_central(i);
		w12 = ext->fn[1]->get_val_neighbor(i);

		w21 = ext->fn[2]->get_val_central(i);
		w22 = ext->fn[2]->get_val_neighbor(i);

		w31 = ext->fn[3]->get_val_central(i);
		w32 = ext->fn[3]->get_val_neighbor(i);

		w_l[0] = w01;
		w_l[1] = w11;
		w_l[2] = w21;
		w_l[3] = w31;

		w_r[0] = w02;
		w_r[1] = w12;
		w_r[2] = w22;
		w_r[3] = w32;

		result += wt[i] * v->val[i] * numerical_flux(2,w_l,w_r,e->nx[i], e->ny[i]);

  }
  return result;
}

double linear_form_interface_3(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w02, w11, w12, w21, w22, w31, w32;
	double w_l[4], w_r[4];
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->get_val_central(i);
		w02 = ext->fn[0]->get_val_neighbor(i);
		
		w11 = ext->fn[1]->get_val_central(i);
		w12 = ext->fn[1]->get_val_neighbor(i);

		w21 = ext->fn[2]->get_val_central(i);
		w22 = ext->fn[2]->get_val_neighbor(i);

		w31 = ext->fn[3]->get_val_central(i);
		w32 = ext->fn[3]->get_val_neighbor(i);

		double p = calc_pressure<double>(w01, w11, w21,w31);

		w_l[0] = w01;
		w_l[1] = w11;
		w_l[2] = w21;
		w_l[3] = w31;

		w_r[0] = w02;
		w_r[1] = w12;
		w_r[2] = w22;
		w_r[3] = w32;

		result += wt[i] * v->val[i] * numerical_flux(3,w_l,w_r,e->nx[i], e->ny[i]);

  }
  return result;
}


///////////////////////////////////////////////////
////////Volumetric linear forms////////////////////
///////////////////////////////////////////////////

double linear_form(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return int_u_v<double,double>(n, wt, ext->fn[0], v) / TAU;
}

///////////////////////////////////////////////////
////////Surface linear forms - IMP///////////////////
///////////////////////////////////////////////////

double linear_form_IMP_0(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w11, w21, w31;
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->val[i];
		
		w11 = ext->fn[1]->val[i];

		w21 = ext->fn[2]->val[i];

		w31 = ext->fn[3]->val[i];

		double a_l = calc_a(w01, w11, w21, w31);
		double a_b = a_l + (R/c_v) * w11 / (2*w01);
		double rho_b = std::pow((a_b*a_b*w01/((1+R/c_v)*calc_pressure<double>(w01,w11,w21,w31))),c_v/R) * w01;
		double p_b = rho_b * a_b * a_b / (R/c_v + 1);
		
		//Ondrej's code.
		double flux[4];
		double alpha = atan2(e->ny[i], e->nx[i]);
		double mat_rot_inv[4][4];
		double flux_local[4];
		flux_local[0] = 0;
		flux_local[1] = p_b;
		flux_local[2] = 0;
		flux_local[3] = 0;
		T_rot(mat_rot_inv, -alpha);
		dot_vector(flux, mat_rot_inv, flux_local);

		result += wt[i] * v->val[i] * flux[0];
  }
  return result;
}

double linear_form_IMP_1(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w11, w21, w31;
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->val[i];
		
		w11 = ext->fn[1]->val[i];

		w21 = ext->fn[2]->val[i];

		w31 = ext->fn[3]->val[i];

		double a_l = calc_a(w01, w11, w21, w31);
		double a_b = a_l + (R/c_v) * w11 / (2*w01);
		double rho_b = std::pow((a_b*a_b*w01/((1+R/c_v)*calc_pressure<double>(w01,w11,w21,w31))),c_v/R) * w01;
		double p_b = rho_b * a_b * a_b / (R/c_v + 1);
		
		//Ondrej's code.
		double flux[4];
		double alpha = atan2(e->ny[i], e->nx[i]);
		double mat_rot_inv[4][4];
		double flux_local[4];
		flux_local[0] = 0;
		flux_local[1] = p_b;
		flux_local[2] = 0;
		flux_local[3] = 0;
		T_rot(mat_rot_inv, -alpha);
		dot_vector(flux, mat_rot_inv, flux_local);

		result += wt[i] * v->val[i] * flux[1];
  }
  return result;
}

double linear_form_IMP_2(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  for (int i = 0; i < n; i++) 
	{

		double w01, w11, w21, w31;

		w01 = ext->fn[0]->val[i];
		
		w11 = ext->fn[1]->val[i];

		w21 = ext->fn[2]->val[i];

		w31 = ext->fn[3]->val[i];

		
		double a_l = calc_a(w01, w11, w21, w31);
		double a_b = a_l + (R/c_v) * w11 / (2*w01);
		double rho_b = std::pow((a_b*a_b*w01/((1+R/c_v)*calc_pressure<double>(w01,w11,w21,w31))),c_v/R) * w01;
		double p_b = rho_b * a_b * a_b / (R/c_v + 1);
		
		//Ondrej's code.
		double flux[4];
		double alpha = atan2(e->ny[i], e->nx[i]);
		double mat_rot_inv[4][4];
		double flux_local[4];
		flux_local[0] = 0;
		flux_local[1] = p_b;
		flux_local[2] = 0;
		flux_local[3] = 0;
		T_rot(mat_rot_inv, -alpha);
		dot_vector(flux, mat_rot_inv, flux_local);
		result += wt[i] * v->val[i] * flux[2];
		
		// Mirroring the state.
		/*
		double wl[4];
		wl[0] = ext->fn[0]->val[i];
		wl[1] = ext->fn[1]->val[i];
		wl[2] = ext->fn[2]->val[i];
		wl[3] = ext->fn[3]->val[i];

		result += wt[i] * v->val[i] * numerical_flux(2,wl,wl,e->nx[i], e->ny[i]);
		*/
  }
  return result;
}

double linear_form_IMP_3(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double w01, w11, w21, w31;
  for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->val[i];
		
		w11 = ext->fn[1]->val[i];

		w21 = ext->fn[2]->val[i];

		w31 = ext->fn[3]->val[i];

		double a_l = calc_a(w01, w11, w21, w31);
		double a_b = a_l + (R/c_v) * w11 / (2*w01);
		double rho_b = std::pow((a_b*a_b*w01/((1+R/c_v)*calc_pressure<double>(w01,w11,w21,w31))),c_v/R) * w01;
		double p_b = rho_b * a_b * a_b / (R/c_v + 1);
		
		//Ondrej's code.
		double flux[4];
		double alpha = atan2(e->ny[i], e->nx[i]);
		double mat_rot_inv[4][4];
		double flux_local[4];
		flux_local[0] = 0;
		flux_local[1] = p_b;
		flux_local[2] = 0;
		flux_local[3] = 0;
		T_rot(mat_rot_inv, -alpha);
		dot_vector(flux, mat_rot_inv, flux_local);

		result += wt[i] * v->val[i] * flux[3];
  }
  return result;
}


///////////////////////////////////////////////////
////////Surface linear forms - IO///////////////////
///////////////////////////////////////////////////

double linear_form_IO_0(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
	double wl[4];
	double wr[4];

  for (int i = 0; i < n; i++) 
	{
		wl[0] = ext->fn[0]->val[i];
		
		wl[1] = ext->fn[1]->val[i];

		wl[2] = ext->fn[2]->val[i];

		wl[3] = ext->fn[3]->val[i];

		if(e->nx[i] < 0)
		{
			double rho_b = bc_density(e->y[i]);
			double u_b = bc_density_vel1(e->y[i]) / bc_density(e->y[i]);
			double v_b = bc_density_vel2(e->y[i]) / bc_density(e->y[i]);

			double a_l = calc_a(wl[0], wl[1], wl[2], wl[3]);
			double a_1 = a_l + (R/c_v) * (wl[1]/wl[0] - u_b);
			
			double rho_1 = std::pow(a_1*a_1*wl[0]/((R/c_v + 1)*calc_pressure(wl[0],wl[1],wl[2],wl[3])), c_v/R) * wl[0];

			double u_1 = u_b;
			double v_1 = wl[2] / wl[0];
			double p_b = rho_1 * a_1 * a_1 / (1+R/c_v);

			double e_1 = calc_energy<double>(rho_1, u_1* rho_1, v_1 * rho_1, p_b);

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * a_l / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/a_l,2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			if(u_b < a_1)
			{
				wr[0] = rho_1;
				wr[1] = u_1 * rho_1;
				wr[2] = v_1 * rho_1;
				wr[3] = e_1;
				result += wt[i] * v->val[i] * numerical_flux(0,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(0,wl,wr,e->nx[i], e->ny[i]);
			}
		}
		else
		{
			double p_b = bc_pressure(e->y[i]);

			double rho_b = wl[0] * std::pow(p_b/calc_pressure(wl[0],wl[1],wl[2],wl[3]),(1/(1+R/c_v)));
			double u_b = (wl[1] / wl[0]) + 2*(c_v/R)*(calc_a<double>(wl[0],wl[1],wl[2],wl[3]) - std::sqrt((1+R/c_v) * p_b / rho_b));
			double v_b = wl[2] / wl[0];
			double e_b = calc_energy<double>(rho_b, u_b*rho_b, v_b*rho_b, p_b);
		

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * calc_a<double>(wl[0],wl[1],wl[2],wl[3]) / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/calc_a<double>(wl[0],wl[1],wl[2],wl[3]),2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			double a_b = calc_a(rho_b, u_b*rho_b, v_b*rho_b, e_b);
			if(u_b < a_b)
			{
				wr[0] = rho_b;
				wr[1] = u_b * rho_b;
				wr[2] = v_b * rho_b;
				wr[3] = e_b;
				result += wt[i] * v->val[i] * numerical_flux(0,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(0,wl,wr,e->nx[i], e->ny[i]);
			}
		}
  }
  return result;
}

double linear_form_IO_1(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
 double result = 0;
	double wl[4];
	double wr[4];

  for (int i = 0; i < n; i++) 
	{
		wl[0] = ext->fn[0]->val[i];
		
		wl[1] = ext->fn[1]->val[i];

		wl[2] = ext->fn[2]->val[i];

		wl[3] = ext->fn[3]->val[i];

		if(e->nx[i] < 0)
		{
			double rho_b = bc_density(e->y[i]);
			double u_b = bc_density_vel1(e->y[i]) / bc_density(e->y[i]);
			double v_b = bc_density_vel2(e->y[i]) / bc_density(e->y[i]);

			double a_l = calc_a(wl[0], wl[1], wl[2], wl[3]);
			double a_1 = a_l + (R/c_v) * (wl[1]/wl[0] - u_b);
			
			double rho_1 = std::pow(a_1*a_1*wl[0]/((R/c_v + 1)*calc_pressure(wl[0],wl[1],wl[2],wl[3])), c_v/R) * wl[0];

			double u_1 = u_b;
			double v_1 = wl[2] / wl[0];
			double p_b = rho_1 * a_1 * a_1 / (1+R/c_v);

			double e_1 = calc_energy<double>(rho_1, u_1* rho_1, v_1 * rho_1, p_b);

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * a_l / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/a_l,2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			if(u_b < a_1)
			{
				wr[0] = rho_1;
				wr[1] = u_1 * rho_1;
				wr[2] = v_1 * rho_1;
				wr[3] = e_1;
				result += wt[i] * v->val[i] * numerical_flux(1,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(1,wl,wr,e->nx[i], e->ny[i]);
			}
		}
		else
		{
			double p_b = bc_pressure(e->y[i]);

			double rho_b = wl[0] * std::pow(p_b/calc_pressure(wl[0],wl[1],wl[2],wl[3]),(1/(1+R/c_v)));
			double u_b = (wl[1] / wl[0]) + 2*(c_v/R)*(calc_a<double>(wl[0],wl[1],wl[2],wl[3]) - std::sqrt((1+R/c_v) * p_b / rho_b));
			double v_b = wl[2] / wl[0];
			double e_b = calc_energy<double>(rho_b, u_b*rho_b, v_b*rho_b, p_b);
		

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * calc_a<double>(wl[0],wl[1],wl[2],wl[3]) / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/calc_a<double>(wl[0],wl[1],wl[2],wl[3]),2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			double a_b = calc_a(rho_b, u_b*rho_b, v_b*rho_b, e_b);
			if(u_b < a_b)
			{
				wr[0] = rho_b;
				wr[1] = u_b * rho_b;
				wr[2] = v_b * rho_b;
				wr[3] = e_b;
				result += wt[i] * v->val[i] * numerical_flux(1,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(1,wl,wr,e->nx[i], e->ny[i]);
			}
		}
  }
  return result;
}

double linear_form_IO_2(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
 double result = 0;
	double wl[4];
	double wr[4];

  for (int i = 0; i < n; i++) 
	{
		wl[0] = ext->fn[0]->val[i];
		
		wl[1] = ext->fn[1]->val[i];

		wl[2] = ext->fn[2]->val[i];

		wl[3] = ext->fn[3]->val[i];

		if(e->nx[i] < 0)
		{
			double rho_b = bc_density(e->y[i]);
			double u_b = bc_density_vel1(e->y[i]) / bc_density(e->y[i]);
			double v_b = bc_density_vel2(e->y[i]) / bc_density(e->y[i]);

			double a_l = calc_a(wl[0], wl[1], wl[2], wl[3]);
			double a_1 = a_l + (R/c_v) * (wl[1]/wl[0] - u_b);
			
			double rho_1 = std::pow(a_1*a_1*wl[0]/((R/c_v + 1)*calc_pressure(wl[0],wl[1],wl[2],wl[3])), c_v/R) * wl[0];

			double u_1 = u_b;
			double v_1 = wl[2] / wl[0];
			double p_b = rho_1 * a_1 * a_1 / (1+R/c_v);

			double e_1 = calc_energy<double>(rho_1, u_1* rho_1, v_1 * rho_1, p_b);

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * a_l / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/a_l,2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			if(u_b < a_1)
			{
				wr[0] = rho_1;
				wr[1] = u_1 * rho_1;
				wr[2] = v_1 * rho_1;
				wr[3] = e_1;
				result += wt[i] * v->val[i] * numerical_flux(2,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(2,wl,wr,e->nx[i], e->ny[i]);
			}
		}
		else
		{
			double p_b = bc_pressure(e->y[i]);

			double rho_b = wl[0] * std::pow(p_b/calc_pressure(wl[0],wl[1],wl[2],wl[3]),(1/(1+R/c_v)));
			double u_b = (wl[1] / wl[0]) + 2*(c_v/R)*(calc_a<double>(wl[0],wl[1],wl[2],wl[3]) - std::sqrt((1+R/c_v) * p_b / rho_b));
			double v_b = wl[2] / wl[0];
			double e_b = calc_energy<double>(rho_b, u_b*rho_b, v_b*rho_b, p_b);
		

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * calc_a<double>(wl[0],wl[1],wl[2],wl[3]) / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/calc_a<double>(wl[0],wl[1],wl[2],wl[3]),2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			double a_b = calc_a(rho_b, u_b*rho_b, v_b*rho_b, e_b);
			if(u_b < a_b)
			{
				wr[0] = rho_b;
				wr[1] = u_b * rho_b;
				wr[2] = v_b * rho_b;
				wr[3] = e_b;
				result += wt[i] * v->val[i] * numerical_flux(2,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(2,wl,wr,e->nx[i], e->ny[i]);
			}
		}
  }
  return result;
}

double linear_form_IO_3(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
 double result = 0;
	double wl[4];
	double wr[4];

  for (int i = 0; i < n; i++) 
	{
		wl[0] = ext->fn[0]->val[i];
		
		wl[1] = ext->fn[1]->val[i];

		wl[2] = ext->fn[2]->val[i];

		wl[3] = ext->fn[3]->val[i];

		if(e->nx[i] < 0)
		{
			double rho_b = bc_density(e->y[i]);
			double u_b = bc_density_vel1(e->y[i]) / bc_density(e->y[i]);
			double v_b = bc_density_vel2(e->y[i]) / bc_density(e->y[i]);

			double a_l = calc_a(wl[0], wl[1], wl[2], wl[3]);
			double a_1 = a_l + (R/c_v) * (wl[1]/wl[0] - u_b);
			
			double rho_1 = std::pow(a_1*a_1*wl[0]/((R/c_v + 1)*calc_pressure(wl[0],wl[1],wl[2],wl[3])), c_v/R) * wl[0];

			double u_1 = u_b;
			double v_1 = wl[2] / wl[0];
			double p_b = rho_1 * a_1 * a_1 / (1+R/c_v);

			double e_1 = calc_energy<double>(rho_1, u_1* rho_1, v_1 * rho_1, p_b);

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * a_l / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/a_l,2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			if(u_b < a_1)
			{
				wr[0] = rho_1;
				wr[1] = u_1 * rho_1;
				wr[2] = v_1 * rho_1;
				wr[3] = e_1;
				result += wt[i] * v->val[i] * numerical_flux(3,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(3,wl,wr,e->nx[i], e->ny[i]);
			}
		}
		else
		{
			double p_b = bc_pressure(e->y[i]);

			double rho_b = wl[0] * std::pow(p_b/calc_pressure(wl[0],wl[1],wl[2],wl[3]),(1/(1+R/c_v)));
			double u_b = (wl[1] / wl[0]) + 2*(c_v/R)*(calc_a<double>(wl[0],wl[1],wl[2],wl[3]) - std::sqrt((1+R/c_v) * p_b / rho_b));
			double v_b = wl[2] / wl[0];
			double e_b = calc_energy<double>(rho_b, u_b*rho_b, v_b*rho_b, p_b);
		

			double a_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * calc_a<double>(wl[0],wl[1],wl[2],wl[3]) / (2+R/c_v);
			double rho_l_star = std::pow(a_l_star/calc_a<double>(wl[0],wl[1],wl[2],wl[3]),2*c_v / R) * wl[0];
			double u_l_star = a_l_star;
			double v_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * a_l_star * a_l_star / (1+R/c_v);
			double e_l_star = calc_energy<double>(rho_l_star, u_l_star * rho_l_star, v_l_star * rho_l_star, p_l_star);

			double a_b = calc_a(rho_b, u_b*rho_b, v_b*rho_b, e_b);
			if(u_b < a_b)
			{
				wr[0] = rho_b;
				wr[1] = u_b * rho_b;
				wr[2] = v_b * rho_b;
				wr[3] = e_b;
				result += wt[i] * v->val[i] * numerical_flux(3,wl,wr,e->nx[i], e->ny[i]);
			}
			else
			{
				wr[0] = rho_l_star;
				wr[1] = u_l_star * rho_l_star;
				wr[2] = v_l_star * rho_l_star;
				wr[3] = e_l_star;
				result += wt[i] * v->val[i] * numerical_flux(3,wl,wr,e->nx[i], e->ny[i]);
			}
		}
  }
  return result;
}