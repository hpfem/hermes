#include "util.h"
#include "limits.h"
#include <limits>

// Calculates energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa)
{
  double to_return = pressure/(kappa - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = (kappa - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = std::sqrt(kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy, kappa) / rho);
  if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

void MachNumberFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(kappa * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / values.at(0)[i]);
}

void PressureFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (kappa - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}

void VelocityFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = values.at(1)[i] / values.at(0)[i];
}

void EntropyFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / p_ext)
      / Hermes::pow((values.at(0)[i] / rho_ext), kappa));
}
