#include "numerical_flux.h"

NumericalFlux::NumericalFlux(double kappa) : kappa(kappa)
{
}

void NumericalFlux::Q(double result[4], double state_vector[4], double nx, double ny)
{
  result[0] = state_vector[0];
  double temp_result_1 = nx * state_vector[1] + ny * state_vector[2];
  double temp_result_2 = -ny * state_vector[1] + nx * state_vector[2];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = state_vector[3];
}

void NumericalFlux::Q_inv(double result[4], double state_vector[4], double nx, double ny)
{
  result[0] = state_vector[0];
  double temp_result_1 = nx * state_vector[1] - ny * state_vector[2];
  double temp_result_2 = ny * state_vector[1] + nx * state_vector[2];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = state_vector[3];
}

UpWindNumericalFlux::UpWindNumericalFlux(double kappa) : NumericalFlux(kappa) {};


double UpWindNumericalFlux::numerical_flux(int component, double w_L[4], double w_R[4], double nx, double ny)
{
  // The upwind direction is for precision calculated on both elements and it is averaged.
  double velocity_cdot_n_L = w_L[1] * nx + w_L[2] * ny;
  double velocity_cdot_n_R = w_R[1] * nx + w_R[2] * ny;
  double velocity_cdot_n = 0.5 * (velocity_cdot_n_L + velocity_cdot_n_R);

  if (velocity_cdot_n > 0)
    return w_L[component];
  else
    return w_R[component];
}

double UpWindNumericalFlux::numerical_flux_inlet(int component, double w_L[4], double w_B[4],
        double nx, double ny)
{
  return w_B[component];
}

double UpWindNumericalFlux::numerical_flux_outlet(int component, double w_L[4], double w_B[4], double nx, double ny)
{
  return w_L[component];
}