// Space distribution of fire temperature.
template<typename Real>
Real T_fire_x(Real x) {
  return -1./32 * x*x*x + 3./16 * x*x;
}

// Temporal distribution of fire temperature.
template<typename Real>
Real T_fire_t(Real t) {
  if (0 <= t  &&  t <= 100) return 0;
  if (100 <= t  &&  t <= 600) return 980. / 500 * (t - 100.);
  if (600 <= t  &&  t <= 1800) return 980;
  if (1800 <= t  &&  t <= 3000) return 980 - 980. / 1200 * (t - 1800.);
  return 0.;
}

// Fire temperature as function of x and time.
template<typename Real>
Real T_fire(Real x, Real t) {
  return T_fire_x(x) * T_fire_t(t) + 20;
}

// Thermal conductivity of the material.
double lambda(double x, double y, double solution_value) {
  return 1.0;
}


