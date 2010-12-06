#include "numerical_flux.h"

extern NumericalFlux num_flux;

static void calc_pressure_func(int n, Hermes::Tuple<double*> scalars, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (num_flux.kappa - 1.) * (scalars.at(3)[i] - 
				(scalars.at(1)[i]*scalars.at(1)[i] + scalars.at(2)[i]*scalars.at(2)[i])/(2*scalars.at(0)[i])
				);
};

static void calc_u_func(int n, Hermes::Tuple<double*> scalars, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = scalars.at(1)[i]/scalars.at(0)[i];
};

static void calc_w_func(int n, Hermes::Tuple<double*> scalars, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = scalars.at(2)[i]/scalars.at(0)[i];
};
