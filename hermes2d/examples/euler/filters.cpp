static void calc_pressure_func(int n, Tuple<scalar*> scalars, scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = R/c_v * (scalars.at(3)[i] - (scalars.at(1)[i]*scalars.at(1)[i] + scalars.at(2)[i]*scalars.at(2)[i])/(2*scalars.at(0)[i]));
};

static void calc_u_func(int n, Tuple<scalar*> scalars, scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = scalars.at(1)[i]/scalars.at(0)[i];
};

static void calc_w_func(int n, Tuple<scalar*> scalars, scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = scalars.at(2)[i]/scalars.at(0)[i];
};