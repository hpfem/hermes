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
				double velocity_x_diag = 0;
				if (u < 0)
					velocity_x_diag = u;
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
					return velocity_x_diag;
				else if (i == 1 && j == 2)
					return 0;
				else if (i == 1 && j == 3)
					return 0;

				else if (i == 2 && j == 0)
					return 0;
				else if (i == 2 && j == 1)
					return 0;
				else if (i == 2 && j == 2)
					return velocity_x_diag;
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

			// calculates the inverterho flux, for testing purposes
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

// The following boundary conditions are the prescribed values used on 
// the inlet and outlet part of the boundary.
// Density boundary condition.
double bc_density(double y)
{
	return 1.0;
}

// Density * velocity in the x coordinate boundary condition.
double bc_density_vel_x(double y)
{
	return 200.0;
}

// Density * velocity in the y coordinate boundary condition.
double bc_density_vel_y(double y)
{
	return 0;
}

// Energy boundary condition.
double bc_energy(double y)
{
	return 1E5;
}
// Calculation of the pressure on the boundary.
double bc_pressure(double y)
{
	return R/c_v * (bc_energy(y) - (bc_density_vel_x(y) * bc_density_vel_x(y) + bc_density_vel_y(y) * bc_density_vel_y(y))/(2 * bc_density(y)));
}


// Initial conditions. Initial conditions correspond with the boundary ones, that is why functions
// specifying boundary conditions are used.
double ic_density(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density(y);
}
double ic_density_vel_x(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density_vel_x(y);
}
double ic_density_vel_y(double x, double y, scalar& dx, scalar& dy)
{
	return bc_density_vel_y(y);
}
double ic_energy(double x, double y, scalar& dx, scalar& dy)
{
	return bc_energy(y);
}

// Calculates pressure from other quantities.
template<typename Scalar>
Scalar calc_pressure(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return (R/c_v) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2*rho));
}
// Calculates the local speed of sound from other quantities.
template<typename Scalar>
Scalar calc_sound_speed(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return std::sqrt(kappa * calc_pressure(rho,rho_v_x,rho_v_y,energy) / rho);
}
// Calculates energy from other quantities.
template<typename Scalar>
Scalar calc_energy(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar pressure)
{
	return (c_v/R) * pressure + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
}

///////////////////////////////////////////////////
////////First flux jacobian/////////////////////
///////////////////////////////////////////////////
template<typename Scalar>
Scalar A_1_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 1;
}

template<typename Scalar>
Scalar A_1_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (R / c_v) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / 	(rho * rho));
}

template<typename Scalar>
Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 2 * (rho_v_x / rho) - (R / c_v) * (rho_v_x / rho);
}

template<typename Scalar>
Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (R / c_v) * (rho_v_y / rho);;
}

template<typename Scalar>
Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return (R / c_v);
}

template<typename Scalar>
Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - rho_v_x * rho_v_y / (rho * rho);
}

template<typename Scalar>
Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_y / rho;
}

template<typename Scalar>
Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_x / rho;
}

template<typename Scalar>
Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (rho_v_x * energy) / (rho * rho) - (rho_v_x / (rho * rho)) * (R / c_v) * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_x / rho) * (R / c_v) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho));
}

template<typename Scalar>
Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return (energy / rho) + (1 / rho) * (R / c_v) * ( energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho))) - (R / c_v) * ((rho_v_x * rho_v_x) / (rho * rho));
}

template<typename Scalar>
Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (R / c_v) * (rho_v_x * rho_v_y) / (rho * rho);
}

template<typename Scalar>
Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_x / rho + (R / c_v) * (rho_v_x / rho);
}


///////////////////////////////////////////////////
////////First flux////////////////////////////////
///////////////////////////////////////////////////

template<typename Scalar>
Scalar f_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return A_1_0_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
		+ A_1_0_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
		+ A_1_0_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
		+ A_1_0_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
}

template<typename Scalar>
Scalar f_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return A_1_1_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
		+ A_1_1_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
		+ A_1_1_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
		+ A_1_1_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
}

template<typename Scalar>
Scalar f_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return A_1_2_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
		+ A_1_2_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
		+ A_1_2_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
		+ A_1_2_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
}

template<typename Scalar>
Scalar f_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return A_1_3_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
		+ A_1_3_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
		+ A_1_3_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
		+ A_1_3_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
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
Scalar A_2_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 1;
}

template<typename Scalar>
Scalar A_2_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - rho_v_x * rho_v_y / (rho * rho);
}

template<typename Scalar>
Scalar A_2_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_y / rho;
}

template<typename Scalar>
Scalar A_2_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_x / rho;
}

template<typename Scalar>
Scalar A_2_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 0;
}

template<typename Scalar>
Scalar A_2_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (R / c_v) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / 	(rho * rho));

}

template<typename Scalar>
Scalar A_2_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (R / c_v) * (rho_v_x / rho);
	
}

template<typename Scalar>
Scalar A_2_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return 2 * (rho_v_y / rho) - (R / c_v) * (rho_v_y / rho);
}

template<typename Scalar>
Scalar A_2_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return R / c_v;
}

template<typename Scalar>
Scalar A_2_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (R / c_v) * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) * (R / c_v) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho));
}

template<typename Scalar>
Scalar A_2_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return - (R / c_v) * (rho_v_x * rho_v_y) / (rho * rho);
}

template<typename Scalar>
Scalar A_2_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return (energy / rho) + (1 / rho) * (R / c_v) * ( energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho))) - (R / c_v) * ((rho_v_y * rho_v_y) / (rho * rho));
}

template<typename Scalar>
Scalar A_2_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy)
{
	return rho_v_y / rho + (R / c_v) * (rho_v_y / rho);
}


// Linear forms coming from time discretization.
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0_time(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * u->val[i] * v->val[i];
	return result / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1_time(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * u->val[i] * v->val[i];
	return result / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2_2_time(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * u->val[i] * v->val[i];
	return result / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_3_3_time(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * u->val[i] * v->val[i];
	return result / TAU;
}

// Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices from the previous time step.
// template : linear_form_m_n means that it is the linear form that it is a bilinear function of a basis function from the m-th space and
// the n-th component of the previous time level solution.
// linear_form_m_n_first_flux and linear_form_m_n_second_flux is distinguished for the forms to be truly bilinear in the above sense.
template<typename Real, typename Scalar>
Scalar linear_form_0_1(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_0_2(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_0_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_0_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_1_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_1_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_2_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_2_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_3_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1_3_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_0_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_0_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_1_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_1_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_2_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_2_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_3_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_2_3_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_0_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_0_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[0]->val[i] * A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_1_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_1_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[1]->val[i] * A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_2_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_2_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[2]->val[i] * A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_3_first_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
	return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_3_3_second_flux(int n, double *wt, Func<Real> *ue[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
	Scalar result = 0;
		for (int i = 0; i < n; i++)
		result += wt[i] * ext->fn[3]->val[i] * A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
	return result;
}

// This is a hack, because of the difficult forms that are used, we supply this artificial integration order. 
Ord	 linear_form_order(int n, double *wt, Func<Ord> *ue[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
	return Ord(20);
}


// Linear forms coming from the DG formulation, evaluated on the inner edges of the triangulation. Linear with respect to the test function v.
// Forms use riemann solvers (where the supplied states are from the previous time level).
// The first function does the calculation, the next ones are added to the weakform and they call the first function, asking for
// the correct element of the vector.
double linear_form_interface(int element, int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
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

		result += wt[i] * v->val[i] * numerical_flux(element,w_l,w_r,e->nx[i], e->ny[i]);
	}
	return result;
}



double linear_form_interface_0(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return linear_form_interface(0, n, wt, ue, v, e, ext);
}
double linear_form_interface_1(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return linear_form_interface(1, n, wt, ue, v, e, ext);
}
double linear_form_interface_2(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return linear_form_interface(2, n, wt, ue, v, e, ext);
}
double linear_form_interface_3(int n, double *wt, Func<double> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return linear_form_interface(3, n, wt, ue, v, e, ext);
}
// Volumetric linear forms. Coming from the time discretization.
// One function used for all the components.
double linear_form(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return int_u_v<double,double>(n, wt, ext->fn[0], v) / TAU;
}

// Surface linear forms representing the solid part of the boundary.
// The flux in the local coordinates is (0, p_b, 0, 0) where p_b stands for pressure on the boundary.
// The first function does the calculation, the next ones are added to the weakform and they call the first function, asking for
// the correct element of the vector.
double bdy_flux_solid_wall_comp(int element, int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	double result = 0;
	double w01, w11, w21, w31;
	for (int i = 0; i < n; i++) 
	{
		w01 = ext->fn[0]->val[i];
		
		w11 = ext->fn[1]->val[i];

		w21 = ext->fn[2]->val[i];

		w31 = ext->fn[3]->val[i];

		double sound_speed_l = calc_sound_speed(w01, w11, w21, w31);
		double sound_speed_b = sound_speed_l + (R/c_v) * w11 / (2*w01);
		double rho_b = std::pow((sound_speed_b*sound_speed_b*w01/(kappa*calc_pressure<double>(w01,w11,w21,w31))),c_v/R) * w01;
		double p_b = rho_b * sound_speed_b * sound_speed_b / kappa;
		
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

		// Mirroring the state
		double w_l[4];
		double w_r[4];
		w_l[0] = w_r[0] = w01;
		w_l[1] = w_r[1] = w11;
		w_l[2] = w_r[2] = w21;
		w_l[3] = w_r[3] = w31;

		result += wt[i] * v->val[i] * numerical_flux(element,w_l,w_r,e->nx[i], e->ny[i]);

		//result += wt[i] * v->val[i] * flux[element];
  }
  return result;
}

double bdy_flux_solid_wall_comp_0(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_solid_wall_comp(0, n, wt, ue, v, e, ext);
}
double bdy_flux_solid_wall_comp_1(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_solid_wall_comp(1, n, wt, ue, v, e, ext);
}
double bdy_flux_solid_wall_comp_2(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_solid_wall_comp(2, n, wt, ue, v, e, ext);
}
double bdy_flux_solid_wall_comp_3(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_solid_wall_comp(3, n, wt, ue, v, e, ext);
}
// Surface linear forms representing the inlet/outlet part of the boundary.
// For details, see comments in the following function.
// The first function does the calculation, the next ones are added to the weakform and they call the first function, asking for
// the correct element of the vector.
double bdy_flux_inlet_outlet_comp(int element, int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	double result = 0;
	// Left (inner) state.
	double wl[4];
	// Eulerian flux.
	double flux[4];

	for (int i = 0; i < n; i++) 
	{
		// Left (inner) state from the previous time level solution.
		wl[0] = ext->fn[0]->val[i];
		
		wl[1] = ext->fn[1]->val[i];

		wl[2] = ext->fn[2]->val[i];

		wl[3] = ext->fn[3]->val[i];

		// The inlet part (left part of the computational domain).
		if(e->nx[i] < 0)
		{
			// Boundary state calculation.
			double rho_b = bc_density(e->y[i]);
			double velocity_x_b = bc_density_vel_x(e->y[i]) / bc_density(e->y[i]);
			double velocity_y_b = bc_density_vel_y(e->y[i]) / bc_density(e->y[i]);

			// Sound speed on the left (inner) side of the boundary.
			double sound_speed_l = calc_sound_speed(wl[0], wl[1], wl[2], wl[3]);

			// Intersection state calculation (marked with an underscore1 (_1)).
			double sound_speed_1 = sound_speed_l + (R/c_v) * (wl[1]/wl[0] - velocity_x_b);
			double rho_1 = std::pow(sound_speed_1*sound_speed_1*wl[0]/(kappa*calc_pressure(wl[0],wl[1],wl[2],wl[3])), c_v/R) * wl[0];
			double velocity_x_1 = velocity_x_b;
			double velocity_y_1 = wl[2] / wl[0];

			// Boundary pressure calculated from the intersection state.
			double p_b = rho_1 * sound_speed_1 * sound_speed_1 / kappa;
			// Calculation of the energy component of the intersection state.
			double energy_1 = calc_energy<double>(rho_1, velocity_x_1* rho_1, velocity_y_1 * rho_1, p_b);

			// Calculation of the state for inflow/outlow velocities above the local speed of sound.
			double sound_speed_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * sound_speed_l / (2+R/c_v);
			double rho_l_star = std::pow(sound_speed_l_star/sound_speed_l,2*c_v / R) * wl[0];
			double velocity_x_l_star = sound_speed_l_star;
			double velocity_y_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * sound_speed_l_star * sound_speed_l_star / kappa;
			double energy_l_star = calc_energy<double>(rho_l_star, velocity_x_l_star * rho_l_star, velocity_y_l_star * rho_l_star, p_l_star);

			// Inflow velocity below the local speed of sound (of the intersection state).
			if(velocity_x_b < sound_speed_1)
			{
				//Ondrej's code.
				double alpha = atan2(e->ny[i], e->nx[i]);
				double mat_rot_inv[4][4];
				double flux_local[4];
				double flux_global[4];
				double mat_rot[4][4];
				T_rot(mat_rot, alpha);
				flux_global[0] = rho_1;
				flux_global[1] = velocity_x_1 * rho_1;
				flux_global[2] = velocity_y_1 * rho_1;
				flux_global[3] = energy_1;
				dot_vector(flux_local, mat_rot, flux_global);
				flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				T_rot(mat_rot_inv, -alpha);
				dot_vector(flux, mat_rot_inv, flux);
			}
			// Inflow velocity above the local speed of sound (of the intersection state).
			else
			{
				//Ondrej's code.
				double alpha = atan2(e->ny[i], e->nx[i]);
				double mat_rot_inv[4][4];
				double flux_local[4];
				double flux_global[4];
				double mat_rot[4][4];
				T_rot(mat_rot, alpha);
				flux_global[0] = rho_l_star;
				flux_global[1] = velocity_x_l_star * rho_l_star;
				flux_global[2] = velocity_y_l_star * rho_l_star;
				flux_global[3] = energy_l_star;
				dot_vector(flux_local, mat_rot, flux_global);
				flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				T_rot(mat_rot_inv, -alpha);
				dot_vector(flux, mat_rot_inv, flux);
			}
		}
		// The outlet part (the right part of the boundary).
		else
		{
			// These calculations are the same as above.
			double p_b = bc_pressure(e->y[i]);
			double rho_b = wl[0] * std::pow(p_b/calc_pressure(wl[0],wl[1],wl[2],wl[3]),(1/kappa));
			double velocity_x_b = (wl[1] / wl[0]) + 2*(c_v/R)*(calc_sound_speed<double>(wl[0],wl[1],wl[2],wl[3]) - std::sqrt(kappa * p_b / rho_b));
			double velocity_y_b = wl[2] / wl[0];
			double energy_b = calc_energy<double>(rho_b, velocity_x_b*rho_b, velocity_y_b*rho_b, p_b);

			double sound_speed_l_star = R/(c_v * (2+R/c_v)) * wl[1] / wl[0] + 2 * calc_sound_speed<double>(wl[0],wl[1],wl[2],wl[3]) / (2+R/c_v);
			double rho_l_star = std::pow(sound_speed_l_star/calc_sound_speed<double>(wl[0],wl[1],wl[2],wl[3]),2*c_v / R) * wl[0];
			double velocity_x_l_star = sound_speed_l_star;
			double velocity_y_l_star = wl[2] / wl[0];
			double p_l_star = rho_l_star * sound_speed_l_star * sound_speed_l_star / kappa;
			double energy_l_star = calc_energy<double>(rho_l_star, velocity_x_l_star * rho_l_star, velocity_y_l_star * rho_l_star, p_l_star);

			double sound_speed_b = calc_sound_speed(rho_b, velocity_x_b*rho_b, velocity_y_b*rho_b, energy_b);
			//  Inflow velocity below the local speed of sound (of the intersection state).
			if(velocity_x_b < sound_speed_b)
			{
				//Ondrej's code.
				double alpha = atan2(e->ny[i], e->nx[i]);
				double mat_rot_inv[4][4];
				double flux_local[4];
				double flux_global[4];
				double mat_rot[4][4];
				T_rot(mat_rot, alpha);
				flux_global[0] = rho_b;
				flux_global[1] = velocity_x_b * rho_b;
				flux_global[2] = velocity_y_b * rho_b;
				flux_global[3] = energy_b;
				dot_vector(flux_local, mat_rot, flux_global);
				flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				T_rot(mat_rot_inv, -alpha);
				dot_vector(flux, mat_rot_inv, flux);
			}
			//  Outflow velocity above the local speed of sound (of the intersection state).
			else
			{
				//Ondrej's code.
				double alpha = atan2(e->ny[i], e->nx[i]);
				double mat_rot_inv[4][4];
				double flux_local[4];
				double flux_global[4];
				double mat_rot[4][4];
				T_rot(mat_rot, alpha);
				flux_global[0] = rho_l_star;
				flux_global[1] = velocity_x_l_star * rho_l_star;
				flux_global[2] = velocity_y_l_star * rho_l_star;
				flux_global[3] = energy_l_star;
				dot_vector(flux_local, mat_rot, flux_global);
				flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
				T_rot(mat_rot_inv, -alpha);
				dot_vector(flux, mat_rot_inv, flux);
			}
		}
	result += wt[i] * v->val[i] * flux[element];
  }
  return result;
}



double bdy_flux_inlet_outlet_comp_0(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_inlet_outlet_comp(0, n, wt, ue, v, e, ext);
}
double bdy_flux_inlet_outlet_comp_1(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_inlet_outlet_comp(1, n, wt, ue, v, e, ext);
}
double bdy_flux_inlet_outlet_comp_2(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_inlet_outlet_comp(2, n, wt, ue, v, e, ext);
}
double bdy_flux_inlet_outlet_comp_3(int n, double *wt, Func<scalar> *ue[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
	return bdy_flux_inlet_outlet_comp(3, n, wt, ue, v, e, ext);
}
