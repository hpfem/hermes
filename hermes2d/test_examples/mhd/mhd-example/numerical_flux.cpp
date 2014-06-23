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

void NumericalFlux::f_1(double result[4], double state[4])
{
  result[0] = state[1];
  result[1] = state[1] * state[1] / state[0] + QuantityCalculator::calc_pressure(state[0], state[1], state[2], state[3], kappa);
  result[2] = state[2] * state[1] / state[0];
  result[3] = (state[1] / state[0]) * (state[3] + QuantityCalculator::calc_pressure(state[0], state[1], state[2], state[3], kappa));
}


VijayasundaramNumericalFlux::VijayasundaramNumericalFlux(double kappa) : StegerWarmingNumericalFlux(kappa)
{
}

void VijayasundaramNumericalFlux::numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny)
{
  double result_temp[4];
  double w_mean[4];
  w_mean[0] = w_L[0] + w_R[0];
  w_mean[1] = w_L[1] + w_R[1];
  w_mean[2] = w_L[2] + w_R[2];
  w_mean[3] = w_L[3] + w_R[3];
  P_plus(result_temp, w_mean, w_L, nx, ny);
  P_minus(result, w_mean, w_R, nx, ny);
  for(unsigned int i = 0; i < 4; i++)
    result[i] += result_temp[i];
}

StegerWarmingNumericalFlux::StegerWarmingNumericalFlux(double kappa) : NumericalFlux(kappa) {};


void StegerWarmingNumericalFlux::numerical_flux(double result[4], double w_L[4], double w_R[4],
        double nx, double ny)
{
  double result_temp[4];
  P_plus(result_temp, w_L, w_L, nx, ny);
  P_minus(result, w_R, w_R, nx, ny);
  for(unsigned int i = 0; i < 4; i++)
    result[i] += result_temp[i];
}
  
double StegerWarmingNumericalFlux::numerical_flux_i(int component, double w_L[4], double w_R[4],
          double nx, double ny)
{
  double result[4];
  numerical_flux(result, w_L, w_R, nx, ny);
  return result[component];
}

void StegerWarmingNumericalFlux::P_plus(double* result, double w[4], double param[4],
          double nx, double ny)
{
  Q(q, w, nx, ny);

  // Initialize the matrices.
  double T[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      T[i][j] = 0.0;

  double T_inv[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      T_inv[i][j] = 0.0;

  // Calculate Lambda^+.
  Lambda_plus(result);

  // Calculate the necessary rows / columns of T(T_inv).
  if(result[0] > 0) {
    T_1(T);
    T_inv_1(T_inv);
  }
  if(result[1] > 0) {
    T_2(T);
    T_inv_2(T_inv);
  }
  if(result[2] > 0) {
    T_3(T);
    T_inv_3(T_inv);
  }
  if(result[3] > 0) {
    T_4(T);
    T_inv_4(T_inv);
  }

  // The matrix T * Lambda * T^{-1}
  double diag_inv[4][4];
  double A_1[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      diag_inv[i][j] = result[i] * T_inv[i][j];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++) {
      A_1[i][j] = 0;
      for(unsigned int k = 0; k < 4; k++)
        A_1[i][j] += T[i][k] * diag_inv[k][j];
    }

  // Finale.
  Q(param, param, nx, ny);
  for(unsigned int i = 0; i < 4; i++) {
    result[i] = 0;
    for(unsigned int j = 0; j < 4; j++)
      result[i] +=A_1[i][j] * param[j];
  }
  Q_inv(result, result, nx, ny);
}

void StegerWarmingNumericalFlux::P_minus(double* result, double w[4], double param[4],
          double nx, double ny)
{
  Q(q, w, nx, ny);

  // Initialize the matrices.
  double T[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      T[i][j] = 0.0;

  double T_inv[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      T_inv[i][j] = 0.0;

  // Calculate Lambda^-.
  Lambda_minus(result);

  // Calculate the necessary rows / columns of T(T_inv).
  if(result[0] < 0) {
    T_1(T);
    T_inv_1(T_inv);
  }
  if(result[1] < 0) {
    T_2(T);
    T_inv_2(T_inv);
  }
  if(result[2] < 0) {
    T_3(T);
    T_inv_3(T_inv);
  }
  if(result[3] < 0) {
    T_4(T);
    T_inv_4(T_inv);
  }


  // The matrix T * Lambda * T^{-1}
  double diag_inv[4][4];
  double A_1[4][4];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      diag_inv[i][j] = result[i] * T_inv[i][j];
  for(unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++) {
      A_1[i][j] = 0;
      for(unsigned int k = 0; k < 4; k++)
        A_1[i][j] += T[i][k] * diag_inv[k][j];
    }

  // Finale.
  Q(param, param, nx, ny);
  for(unsigned int i = 0; i < 4; i++) {
    result[i] = 0;
    for(unsigned int j = 0; j < 4; j++)
      result[i] +=A_1[i][j] * param[j];
  }
  Q_inv(result, result, nx, ny);
}

void StegerWarmingNumericalFlux::Lambda_plus(double result[4])
{
  a = QuantityCalculator::calc_sound_speed(q[0], q[1], q[2], q[3], kappa);
  u = q[1] / q[0];
  v = q[2] / q[0];
  V = u*u + v*v;
  result[0] = u - a < 0 ? 0 : u - a;
  result[1] = u < 0 ? 0 : u;
  result[2] = u < 0 ? 0 : u;
  result[3] = u + a < 0 ? 0 : u + a;
}

void StegerWarmingNumericalFlux::Lambda_minus(double result[4])
{
  a = QuantityCalculator::calc_sound_speed(q[0], q[1], q[2], q[3], kappa);
  u = q[1] / q[0];
  v = q[2] / q[0];
  V = u*u + v*v;
  result[0] = u - a < 0 ? u - a : 0;
  result[1] = u < 0 ? u : 0;
  result[2] = u < 0 ? u : 0;
  result[3] = u + a < 0 ? u + a : 0;
}

void StegerWarmingNumericalFlux::Lambda(double result[4])
{
  a = QuantityCalculator::calc_sound_speed(q[0], q[1], q[2], q[3], kappa);
  u = q[1] / q[0];
  v = q[2] / q[0];
  V = u*u + v*v;
  result[0] = u - a ;
  result[1] = u;
  result[2] = u;
  result[3] = u + a;
}

void StegerWarmingNumericalFlux::T_1(double result[4][4])
{
  result[0][0] = 1.0;
  result[1][0] = u - a;
  result[2][0] = v;
  result[3][0] = (V / 2.0) + (a*a / (kappa - 1.0)) - (u * a);
}
void StegerWarmingNumericalFlux::T_2(double result[4][4])
{
  result[0][1] = 1.0;
  result[1][1] = u;
  result[2][1] = v;
  result[3][1] = V / 2.0;
}
void StegerWarmingNumericalFlux::T_3(double result[4][4])
{
  result[0][2] = 1.0;
  result[1][2] = u;
  result[2][2] = v - a;
  result[3][2] = (V / 2.0)  - v * a;
}
void StegerWarmingNumericalFlux::T_4(double result[4][4])
{
  result[0][3] = 1.0;
  result[1][3] = u + a;
  result[2][3] = v;
  result[3][3] = (V / 2.0) + (a * a / (kappa - 1.0)) + (u * a);
}

void StegerWarmingNumericalFlux::T_inv_1(double result[4][4])
{
  result[0][0] = (1.0 / (a * a)) * (0.5 * (((kappa - 1) * V / 2.0) + u * a));
  result[0][1] = (1.0 / (a * a)) * (- (a + u * (kappa - 1.0)) / 2.0);
  result[0][2] = (1.0 / (a * a)) * (- (v * (kappa - 1.0)) / 2.0);
  result[0][3] = (1.0 / (a * a)) * (kappa - 1.0) / 2.0;
}
void StegerWarmingNumericalFlux::T_inv_2(double result[4][4])
{
  result[1][0] = (1.0 / (a * a)) * (a * a - v * a - (kappa - 1.0) * (V / 2.0));
  result[1][1] = (1.0 / (a * a)) * u * (kappa - 1.0);
  result[1][2] = (1.0 / (a * a)) * (a + v * (kappa - 1.0));
  result[1][3] = (1.0 / (a * a)) * (1.0 - kappa);
}
void StegerWarmingNumericalFlux::T_inv_3(double result[4][4])
{
  result[2][0] = (1.0 / (a * a)) * v * a;
  result[2][1] = (1.0 / (a * a)) * 0.0;
  result[2][2] = (1.0 / (a * a)) * (-a);
  result[2][3] = (1.0 / (a * a)) * 0.0;
}
void StegerWarmingNumericalFlux::T_inv_4(double result[4][4])
{
  result[3][0] = (1.0 / (a * a)) * (0.5 * (((kappa - 1.0) * V / 2.0) - u * a));
  result[3][1] = (1.0 / (a * a)) * (a - u * (kappa - 1.0)) / 2.0;
  result[3][2] = (1.0 / (a * a)) * (- (v * (kappa - 1.0)) / 2.0);
  result[3][3] = (1.0 / (a * a)) * (kappa - 1.0) / 2.0;
}

void StegerWarmingNumericalFlux::numerical_flux_solid_wall(double result[4], double w_L[4], double nx, double ny)
{
  Q(q_L, w_L, nx, ny);
  a_B = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa) + ((kappa - 1) * q_L[1] / (2 * q_L[0]));
  double rho_B = std::pow(a_B * a_B * q_L[0] / (kappa * QuantityCalculator::calc_pressure(q_L[0], q_L[1], q_L[2], q_L[3], kappa)), (1 / (kappa - 1))) * q_L[0];
  q_R[0] = 0;
  q_R[1] = rho_B * a_B * a_B / kappa;
  q_R[2] = 0;
  q_R[3] = 0;
  Q_inv(result, q_R, nx, ny);
}
  
double StegerWarmingNumericalFlux::numerical_flux_solid_wall_i(int component, double w_L[4], double nx, double ny)
{
  double result[4];
  numerical_flux_solid_wall(result, w_L, nx, ny);
  return result[component];
}

void StegerWarmingNumericalFlux::numerical_flux_inlet(double result[4], double w_L[4], double w_B[4],
        double nx, double ny)
{
  // At the beginning, rotate the states into the local coordinate system and store the left and right state
  // so we do not have to pass it around.
  Q(q_L, w_L, nx, ny);
  Q(q_B, w_B, nx, ny);

  // Speeds of sound.
  a_L = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa);
  a_B = QuantityCalculator::calc_sound_speed(q_B[0], q_B[1], q_B[2], q_B[3], kappa);

  if(q_L[1] / q_L[0] > a_L) {// Supersonic inlet - everything is prescribed.
    f_1(result, q_B);
    Q_inv(result, result, nx, ny);
    return;
  }
  else {// Subsonic inlet - only rho_b, v_x, v_y are prescribed, pressure is calculated as follows. The pressure is prescribed always so that one can know if the
    // inlet is subsonic or supersonic.
    double a_1 = a_L + ((kappa - 1) / 2) * (q_L[1] / q_L[0] - q_B[1] / q_B[0]);
    q_1[0] = std::pow(a_1 * a_1 * q_L[0] / (kappa * QuantityCalculator::calc_pressure(q_L[0], q_L[1], q_L[2], q_L[3], kappa)), 1 / (kappa - 1)) * q_L[0];
    q_1[1] = q_1[0] * q_B[1] / q_B[0];
    q_1[2] = q_1[0] * q_L[2] / q_L[0];
    q_1[3] = QuantityCalculator::calc_energy(q_1[0], q_1[1], q_1[2], a_1 * a_1 * q_1[0] / kappa, kappa);
    if(q_B[1] / q_B[0] < 0)
      if(q_B[1] / q_B[0] < a_1) {
        f_1(result, q_B);
        Q_inv(result, result, nx, ny);
        return;
      }
      else {
        double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
        q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
        q_L_star[1] = a_l_star;
        q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
        q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
        double first1[4];
        double second1[4];
        double third1[4];
        f_1(first1, q_B);
        f_1(second1, q_L_star);
        f_1(third1, q_1);
        for(unsigned int i = 0; i < 4; i++)
          result[i] = first1[i] + second1[i] - third1[i];
        Q_inv(result, result, nx, ny);
        return;
      }
    else
      if(q_B[1] / q_B[0] < a_1) {
        f_1(result, q_1);
        Q_inv(result, result, nx, ny);
        return;
      }
      else {
        double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
        q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
        q_L_star[1] = a_l_star;
        q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
        q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
        f_1(result, q_L_star);
        Q_inv(result, result, nx, ny);
        return;
      }
  }
}
  
double StegerWarmingNumericalFlux::numerical_flux_inlet_i(int component, double w_L[4], double w_B[4],
        double nx, double ny)
{
  double result[4];
  numerical_flux_inlet(result, w_L, w_B, nx, ny);
  return result[component];
}

void StegerWarmingNumericalFlux::numerical_flux_outlet(double result[4], double w_L[4], double pressure, double nx, double ny)
{
  // At the beginning, rotate the states into the local coordinate system and store the left and right state
  // so we do not have to pass it around.
  Q(q_L, w_L, nx, ny);

  double a_L = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa);

  if(q_L[1] / q_L[0] > a_L) {// Supersonic inlet - everything is prescribed.
    f_1(result, q_L);
    Q_inv(result, result, nx, ny);
    return;
  }
  else {
    this->q_B[0] = q_L[0] * std::pow(pressure / QuantityCalculator::calc_pressure(this->q_L[0], this->q_L[1], this->q_L[2], this->q_L[3], kappa), 1 / kappa);
    this->q_B[1] = this->q_B[0] * (q_L[1] / q_L[0] + (2 / (kappa - 1)) * (a_L - std::sqrt(kappa * pressure / q_B[0])));
    this->q_B[2] = this->q_B[0] * this->q_L[2] / this->q_L[0];
    this->q_B[3] = QuantityCalculator::calc_energy(this->q_B[0], this->q_B[1], this->q_B[2], pressure, kappa);
    if(q_B[1] / q_B[0] < QuantityCalculator::calc_sound_speed(this->q_B[0], this->q_B[1], this->q_B[2], this->q_B[3], kappa)) {
      f_1(result, q_B);
      Q_inv(result, result, nx, ny);
      return;
    }
    else {
      double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
      q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
      q_L_star[1] = a_l_star;
      q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
      q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
      f_1(result, q_L_star);
      Q_inv(result, result, nx, ny);
      return;
    }
  }
}
  
double StegerWarmingNumericalFlux::numerical_flux_outlet_i(int component, double w_L[4], double pressure, double nx, double ny)
{
  double result[4];
  numerical_flux_outlet(result, w_L, pressure, nx, ny);
  return result[component];
}

double* StegerWarmingNumericalFlux::get_q()
{
  return q;
}



OsherSolomonNumericalFlux::OsherSolomonNumericalFlux(double kappa) : NumericalFlux(kappa)
{
}

void OsherSolomonNumericalFlux::numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny)
{

  // At the beginning, rotate the states into the local coordinate system and store the left and right state
  // so we do not have to pass it around.
  Q(q_L, w_L, nx, ny);
  Q(q_R, w_R, nx, ny);

  // Decide what we have to calculate.
  // Speeds of sound.
  a_L = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa);
  a_R = QuantityCalculator::calc_sound_speed(q_R[0], q_R[1], q_R[2], q_R[3], kappa);
  
  // Check that we can use the following.
  double right_hand_side = 0;
  if((q_L[2] / q_L[0]) - (q_R[2] / q_R[0]) > 0)
    right_hand_side = (q_L[2] / q_L[0] - q_R[2] / q_R[0] > 0) / 2;
  if(a_L + a_R + ((kappa - 1) * (q_L[1] / q_L[0] - q_R[1] / q_R[0]) / 2) <= right_hand_side)
    throw Hermes::Exceptions::Exception("Osher-Solomon numerical flux is not possible to construct according to the table.");

  // Utility numbers.
  this->z_L = (0.5 * (kappa - 1) * q_L[1] / q_L[0]) + a_L;
  this->z_R = (0.5 * (kappa - 1) * q_R[1] / q_R[0]) - a_R;
  this->s_L = QuantityCalculator::calc_pressure(q_L[0], q_L[1], q_L[2], q_L[3], kappa) / std::pow(q_L[0], kappa);
  this->s_R = QuantityCalculator::calc_pressure(q_R[0], q_R[1], q_R[2], q_R[3], kappa) / std::pow(q_R[0], kappa);
  this->alpha = std::pow(s_R / s_L, 1 / (2 * kappa));

  // We always need to calculate q_1, a_1, a_3, as we are going to decide what to return based on this.
  calculate_q_1_a_1_a_3();

  // First column in table 3.4.1 on the page 233 in Feist (2003).
  if(q_R[1] / q_R[0] >= - a_R && q_L[1] / q_L[0] <= a_L) {
    // First row.
    if(a_1 <= q_1[1] / q_1[0]) {
      calculate_q_L_star();
      f_1(result, q_L_star);
      Q_inv(result, result, nx, ny);
      return;
    }
    // Second row.
    if(0 < q_1[1] / q_1[0] && q_1[1] / q_1[0] < a_1) {
      f_1(result, q_1);
      Q_inv(result, result, nx, ny);
      return;
    }
    // Third row.
    if(-a_3 <= q_1[1] / q_1[0] && q_1[1] / q_1[0] <= 0) {
      calculate_q_3();
      f_1(result, q_3);
      Q_inv(result, result, nx, ny);
      return;
    }
    // Fourth row.
    if(q_1[1] / q_1[0] < -a_3) {
      calculate_q_R_star();
      f_1(result, q_R_star);
      Q_inv(result, result, nx, ny);
      return;
    }
  }

  // Second column in table 3.4.1 on the page 233 in Feist (2003).
  if(q_R[1] / q_R[0] >= - a_R && q_L[1] / q_L[0] > a_L) {
    // First row.
    if(a_1 <= q_1[1] / q_1[0]) {
      f_1(result, q_L);
      Q_inv(result, result, nx, ny);
      return;
    }
    // Second row.
    if(0 < q_1[1] / q_1[0] && q_1[1] / q_1[0] < a_1) {
      calculate_q_L_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_L);
      f_1(second1, q_L_star);
      f_1(third1, q_1);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Third row.
    if(-a_3 <= q_1[1] / q_1[0] && q_1[1] / q_1[0] <= 0) {
      calculate_q_L_star();
      calculate_q_3();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_L);
      f_1(second1, q_L_star);
      f_1(third1, q_3);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Fourth row.
    if(q_1[1] / q_1[0] < -a_3) {
      calculate_q_L_star();
      calculate_q_R_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_L);
      f_1(second1, q_L_star);
      f_1(third1, q_R_star);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
  }

  // Third column in table 3.4.1 on the page 233 in Feist (2003).
  if(q_R[1] / q_R[0] < - a_R && q_L[1] / q_L[0] <= a_L) {
    // First row.
    if(a_1 <= q_1[1] / q_1[0]) {
      calculate_q_R_star();
      calculate_q_L_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_R);
      f_1(second1, q_R_star);
      f_1(third1, q_L_star);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Second row.
    if(0 < q_1[1] / q_1[0] && q_1[1] / q_1[0] < a_1) {
      calculate_q_R_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_R);
      f_1(second1, q_R_star);
      f_1(third1, q_1);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Third row.
    if(-a_3 <= q_1[1] / q_1[0] && q_1[1] / q_1[0] <= 0) {
      calculate_q_R_star();
      calculate_q_3();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_R);
      f_1(second1, q_R_star);
      f_1(third1, q_3);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      return;
    }
    // Fourth row.
    if(q_1[1] / q_1[0] < -a_3) {
      f_1(result, q_R);
      Q_inv(result, result, nx, ny);
      return;
    }
  }

  // Fourth column in table 3.4.1 on the page 233 in Feist (2003).
  if(q_R[1] / q_R[0] < - a_R && q_L[1] / q_L[0] > a_L) {
    // First row.
    if(a_1 <= q_1[1] / q_1[0]) {
      calculate_q_R_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_L);
      f_1(second1, q_R_star);
      f_1(third1, q_R);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Second row.
    if(0 < q_1[1] / q_1[0] && q_1[1] / q_1[0] < a_1) {
      calculate_q_R_star();
      calculate_q_L_star();
      double first1[4];
      double second1[4];
      double third1[4];
      double fourth1[4];
      double fifth1[4];
      f_1(first1, q_L);
      f_1(second1, q_R_star);
      f_1(third1, q_R);
      f_1(fourth1, q_L_star);
      f_1(fifth1, q_1);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i] - fourth1[i] + fifth1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Third row.
    if(-a_3 <= q_1[1] / q_1[0] && q_1[1] / q_1[0] <= 0) {
      calculate_q_R_star();
      calculate_q_L_star();
      calculate_q_3();
      double first1[4];
      double second1[4];
      double third1[4];
      double fourth1[4];
      double fifth1[4];
      f_1(first1, q_L);
      f_1(second1, q_R_star);
      f_1(third1, q_R);
      f_1(fourth1, q_L_star);
      f_1(fifth1, q_3);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] - second1[i] + third1[i] - fourth1[i] + fifth1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
    // Fourth row.
    if(q_1[1] / q_1[0] < -a_3) {
      calculate_q_L_star();
      double first1[4];
      double second1[4];
      double third1[4];
      f_1(first1, q_L);
      f_1(second1, q_R);
      f_1(third1, q_L_star);
      for(unsigned int i = 0; i < 4; i++)
        result[i] = first1[i] + second1[i] - third1[i];
      Q_inv(result, result, nx, ny);
      return;
    }
  }
}

void OsherSolomonNumericalFlux::numerical_flux_solid_wall(double result[4], double w_L[4], double nx, double ny)
{
  Q(q_L, w_L, nx, ny);
  a_B = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa) + ((kappa - 1) * q_L[1] / (2 * q_L[0]));
  double rho_B = std::pow(a_B * a_B * q_L[0] / (kappa * QuantityCalculator::calc_pressure(q_L[0], q_L[1], q_L[2], q_L[3], kappa)), (1 / (kappa - 1))) * q_L[0];
  q_R[0] = 0;
  q_R[1] = rho_B * a_B * a_B / kappa;
  q_R[2] = 0;
  q_R[3] = 0;
  Q_inv(result, q_R, nx, ny);
}
  
double OsherSolomonNumericalFlux::numerical_flux_solid_wall_i(int component, double w_L[4], double nx, double ny)
{
  double result[4];
  numerical_flux_solid_wall(result, w_L, nx, ny);
  return result[component];
}

void OsherSolomonNumericalFlux::numerical_flux_inlet(double result[4], double w_L[4], double w_B[4],
        double nx, double ny)
{
  // At the beginning, rotate the states into the local coordinate system and store the left and right state
  // so we do not have to pass it around.
  Q(q_L, w_L, nx, ny);
  Q(q_B, w_B, nx, ny);

  // Speeds of sound.
  a_L = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa);
  a_B = QuantityCalculator::calc_sound_speed(q_B[0], q_B[1], q_B[2], q_B[3], kappa);

  if(q_L[1] / q_L[0] > a_L) {// Supersonic inlet - everything is prescribed.
    f_1(result, q_B);
    Q_inv(result, result, nx, ny);
    return;
  }
  else {// Subsonic inlet - only rho_b, v_x, v_y are prescribed, pressure is calculated as follows. The pressure is prescribed always so that one can know if the
    // inlet is subsonic or supersonic.
    double a_1 = a_L + ((kappa - 1) / 2) * (q_L[1] / q_L[0] - q_B[1] / q_B[0]);
    q_1[0] = std::pow(a_1 * a_1 * q_L[0] / (kappa * QuantityCalculator::calc_pressure(q_L[0], q_L[1], q_L[2], q_L[3], kappa)), 1 / (kappa - 1)) * q_L[0];
    q_1[1] = q_1[0] * q_B[1] / q_B[0];
    q_1[2] = q_1[0] * q_L[2] / q_L[0];
    q_1[3] = QuantityCalculator::calc_energy(q_1[0], q_1[1], q_1[2], a_1 * a_1 * q_1[0] / kappa, kappa);
    if(q_B[1] / q_B[0] < 0)
      if(q_B[1] / q_B[0] < a_1) {
        f_1(result, q_B);
        Q_inv(result, result, nx, ny);
        return;
      }
      else {
        double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
        q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
        q_L_star[1] = a_l_star;
        q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
        q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
        double first1[4];
        double second1[4];
        double third1[4];
        f_1(first1, q_B);
        f_1(second1, q_L_star);
        f_1(third1, q_1);
        for(unsigned int i = 0; i < 4; i++)
          result[i] = first1[i] + second1[i] - third1[i];
        Q_inv(result, result, nx, ny);
        return;
      }
    else
      if(q_B[1] / q_B[0] < a_1) {
        f_1(result, q_1);
        Q_inv(result, result, nx, ny);
        return;
      }
      else {
        double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
        q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
        q_L_star[1] = a_l_star;
        q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
        q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
        f_1(result, q_L_star);
        Q_inv(result, result, nx, ny);
        return;
      }
  }
}
  
double OsherSolomonNumericalFlux::numerical_flux_inlet_i(int component, double w_L[4], double w_B[4],
        double nx, double ny)
{
  double result[4];
  numerical_flux_inlet(result, w_L, w_B, nx, ny);
  return result[component];
}

void OsherSolomonNumericalFlux::numerical_flux_outlet(double result[4], double w_L[4], double pressure, double nx, double ny)
{
  // At the beginning, rotate the states into the local coordinate system and store the left and right state
  // so we do not have to pass it around.
  Q(q_L, w_L, nx, ny);

  double a_L = QuantityCalculator::calc_sound_speed(q_L[0], q_L[1], q_L[2], q_L[3], kappa);

  if(q_L[1] / q_L[0] > a_L) {// Supersonic inlet - everything is prescribed.
    f_1(result, q_L);
    Q_inv(result, result, nx, ny);
    return;
  }
  else {
    this->q_B[0] = q_L[0] * std::pow(pressure / QuantityCalculator::calc_pressure(this->q_L[0], this->q_L[1], this->q_L[2], this->q_L[3], kappa), 1 / kappa);
    this->q_B[1] = this->q_B[0] * (q_L[1] / q_L[0] + (2 / (kappa - 1)) * (a_L - std::sqrt(kappa * pressure / q_B[0])));
    this->q_B[2] = this->q_B[0] * this->q_L[2] / this->q_L[0];
    this->q_B[3] = QuantityCalculator::calc_energy(this->q_B[0], this->q_B[1], this->q_B[2], pressure, kappa);
    if(q_B[1] / q_B[0] < QuantityCalculator::calc_sound_speed(this->q_B[0], this->q_B[1], this->q_B[2], this->q_B[3], kappa)) {
      f_1(result, q_B);
      Q_inv(result, result, nx, ny);
      return;
    }
    else {
      double a_l_star = (((kappa - 1) / (kappa + 1)) * q_L[1] / q_L[0]) + 2 * a_L / (kappa + 1);
      q_L_star[0] = std::pow(a_l_star / a_L, 2 / (kappa - 1)) * q_L[0];
      q_L_star[1] = a_l_star;
      q_L_star[2] = q_L_star[0] * q_L[2] / q_L[0];
      q_L_star[3] = QuantityCalculator::calc_energy(q_L_star[0], q_L_star[1], q_L_star[2], q_L_star[0] * a_l_star * a_l_star / kappa, kappa);
      f_1(result, q_L_star);
      Q_inv(result, result, nx, ny);
      return;
    }
  }
}
  
double OsherSolomonNumericalFlux::numerical_flux_outlet_i(int component, double w_L[4], double pressure, double nx, double ny)
{
  double result[4];
  numerical_flux_outlet(result, w_L, pressure, nx, ny);
  return result[component];
}

void OsherSolomonNumericalFlux::calculate_q_1_a_1_a_3()
{
  this->a_1 = (z_L - z_R) / (1 + alpha);
  this->q_1[0] = std::pow(a_1 / a_L, 2 / (kappa - 1)) * q_L[0];
  this->q_1[1] = this->q_1[0] * 2 * (z_L - a_1) / (kappa - 1);
  this->q_1[2] = this->q_1[0] * this->q_L[2] / this->q_L[0] ;
  this->q_1[3] = QuantityCalculator::calc_energy(this->q_1[0], this->q_1[1], this->q_1[2], a_1 * a_1 * q_1[0] / kappa, kappa);
  this->a_3 = alpha * a_1;
}

void OsherSolomonNumericalFlux::calculate_q_L_star()
{
  this->a_L_star = 2 * z_L / (kappa + 1);
  this->q_L_star[0] = std::pow(a_L_star / a_L, 2 / (kappa -1 )) * q_L[0];
  this->q_L_star[1] = this->q_L_star[0] * a_L_star;
  this->q_L_star[2] = this->q_1[0] * this->q_L[2] / this->q_L[0] ;
  this->q_L_star[3] = QuantityCalculator::calc_energy(this->q_L_star[0], this->q_L_star[1], this->q_L_star[2], a_L_star * a_L_star * q_L_star[0] / kappa, kappa);
}

void OsherSolomonNumericalFlux::calculate_q_3()
{
  // a_3 already calculated.
  this->q_3[0] = q_1[0] / (alpha * alpha);
  this->q_3[1] = this->q_3[0] * this->q_1[1] / this->q_1[0] ;
  this->q_3[2] = this->q_3[0] * this->q_R[2] / this->q_R[0] ;
  this->q_3[3] = QuantityCalculator::calc_energy(this->q_3[0], this->q_3[1], this->q_3[2], a_3 * a_3 * q_3[0] / kappa, kappa);
}

void OsherSolomonNumericalFlux::calculate_q_R_star()
{
  this->a_R_star = - 2 * z_R / (kappa + 1);
  this->q_R_star[0] = std::pow(a_R_star / a_R, 2 / (kappa -1 )) * q_R[0];
  this->q_R_star[1] = this->q_R_star[0] * -a_R_star;
  this->q_R_star[2] = this->q_1[0] * this->q_R[2] / this->q_R[0] ;
  this->q_R_star[3] = QuantityCalculator::calc_energy(this->q_R_star[0], this->q_R_star[1], this->q_R_star[2], a_R_star * a_R_star * q_R_star[0] / kappa, kappa);
}

double OsherSolomonNumericalFlux::numerical_flux_i(int component, double w_L[4], double w_R[4],
          double nx, double ny)
{
  double result[4];
  numerical_flux(result, w_L, w_R, nx, ny);
  return result[component];
}






















/*
double NumericalFlux::f_x(int component, double w0, double w1, double w3, double w4)
{
    if (i == 0)
        return w1;
    else if (i == 1)
        return w1*w1/w0 + (kappa - 1.) * (w4 - (w1*w1+w3*w3)/(2*w0));
    else if (i == 2)
        return w1*w3/w0;
    else if (i == 3)
        return w1/w0 * (w4 + (kappa - 1.) * (w4 - (w1*w1+w3*w3)/(2*w0)));

    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

double NumericalFlux::f_z(int component, double w0, double w1, double w3, double w4)
{
    if (i == 0)
        return w3;
    else if (i == 1)
        return w3*w1/w0;
    else if (i == 2)
        return w3*w3/w0 + (kappa - 1.) * (w4 - (w1*w1+w3*w3)/(2*w0));
    else if (i == 3)
        return w3/w0 * (w4 + (kappa - 1.) * (w4 - (w1*w1+w3*w3)/(2*w0)));

    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

double NumericalFlux::A_x(int component, int j, double w0, double w1, double w3, double w4)
{
    if (i == 0 && j == 0)
        return 0;
    else if (i == 0 && j == 1)
        return 1;
    else if (i == 0 && j == 2)
        return 0;
    else if (i == 0 && j == 3)
        return 0;

    else if (i == 1 && j == 0)
        return -w1*w1/(w0*w0) + (kappa - 1.) * (w1*w1 + w3*w3)/(2 * w0*w0);
    else if (i == 1 && j == 1)
        return 2*w1/w0 - (kappa - 1.) * w1 / w0;
    else if (i == 1 && j == 2)
        return - (kappa - 1.) * w3 / w0;
    else if (i == 1 && j == 3)
        return kappa - 1.;

    else if (i == 2 && j == 0)
        return -w1*w3/(w0*w0);
    else if (i == 2 && j == 1)
        return w3/w0;
    else if (i == 2 && j == 2)
        return w1/w0;
    else if (i == 2 && j == 3)
        return 0;

    else if (i == 3 && j == 0)
        return -w1*w4/(w0*w0) - w1/(w0*w0) * (kappa - 1.)
            * (w4 - (w1*w1+w3*w3)/(2*w0)) + w1/w0 * (kappa - 1.)
            * (w1*w1+w3*w3)/(2*w0*w0);
        // or equivalently:
        //return w1/w0 * ((kappa - 1.) * (w1*w1+w3*w3)/(w0*w0) - ((kappa - 1.) + 1) * w4/w0);
    else if (i == 3 && j == 1)
        return w4/w0 + 1/w0 * (kappa - 1.)
            * (w4 - (w1*w1+w3*w3)/(2*w0)) - (kappa - 1.)
            * w1*w1/(w0*w0);
    else if (i == 3 && j == 2)
        return - (kappa - 1.) * w1*w3/(w0*w0);
    else if (i == 3 && j == 3)
        return w1/w0 + (kappa - 1.) * w1/w0;

    printf("i=%d, j=%d;\n", i, j);
    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

double NumericalFlux::A_z(int component, int j, double w0, double w1, double w3, double w4)
{
    if (i == 0 && j == 0)
        return 0;
    else if (i == 0 && j == 1)
        return 0;
    else if (i == 0 && j == 2)
        return 1;
    else if (i == 0 && j == 3)
        return 0;

    else if (i == 1 && j == 0)
        return -w3*w1/(w0*w0);
    else if (i == 1 && j == 1)
        return w3/w0;
    else if (i == 1 && j == 2)
        return w1/w0;
    else if (i == 1 && j == 3)
        return 0;

    else if (i == 2 && j == 0)
        return -w3*w3/(w0*w0) + (kappa - 1.) * (w1*w1 + w3*w3)/(2*w0*w0);
    else if (i == 2 && j == 1)
        return - (kappa - 1.) * w1 / w0;
    else if (i == 2 && j == 2)
        return 2*w3/w0 - (kappa - 1.) * w3 / w0;
    else if (i == 2 && j == 3)
        return (kappa - 1.);

    else if (i == 3 && j == 0)
        return -w3*w4/(w0*w0) - w3/(w0*w0) * (kappa - 1.)
            * (w4 - (w1*w1+w3*w3)/(2*w0)) + w3/w0 * (kappa - 1.)
            * (w1*w1+w3*w3)/(2*w0*w0);
        // or equivalently:
        //return w1/w0 * ((kappa - 1.) * (w1*w1+w3*w3)/(w0*w0) - ((kappa - 1.) + 1) * w4/w0);
    else if (i == 3 && j == 1)
        return - (kappa - 1.) * w3*w1/(w0*w0);
    else if (i == 3 && j == 2)
        return w4/w0 + 1/w0 * (kappa - 1.)
            * (w4 - (w1*w1+w3*w3)/(2*w0)) - (kappa - 1.)
            * w3*w3/(w0*w0);
    else if (i == 3 && j == 3)
        return w3/w0 + (kappa - 1.) * w3/w0;

    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

double NumericalFlux::matrix_R(int component, int j, double w0, double w1, double w3, double w4)
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
    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

double NumericalFlux::matrix_R_inv(int component, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double p = (kappa-1)*(E - rho*v2/2);
    double c = sqrt(kappa*p/rho);
    double result = 0;
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
        throw Hermes::Exceptions::Exception("Invalid index.");
    }
    return result/(c*c);
}

double NumericalFlux::matrix_D_minus(int component, int j, double w0, double w1, double w3, double w4)
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
    throw Hermes::Exceptions::Exception("Invalid index.");
    return 0.0;
}

// multiplies two matrices
void NumericalFlux::dot(double result[4][4], double A[4][4], double B[4][4])
{
    for (int component=0; i < 4; i++)
        for (int j=0; j < 4; j++) {
            double sum=0;
            for (int k=0; k < 4; k++)
                sum += A[i][k] * B[k][j];
            result[i][j] = sum;
        }
}

// multiplies a matrix and a vector
void NumericalFlux::dot_vector(double result[4], double A[4][4], double B[4])
{
    for (int component=0; i < 4; i++) {
        double sum=0;
        for (int k=0; k < 4; k++)
            sum += A[i][k] * B[k];
        result[i] = sum;
    }
}

// XXX: this matrix should take the normals directly, e.g.
// [cos, sin]
// [-sin, cos]
// becomes
// [nx, ny]
// [-ny, nx]
void NumericalFlux::T_rot(double result[4][4], double beta)
{
    for (int component=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            result[i][j] = 0;
    result[0][0] = 1;
    result[1][1] = cos(beta);
    result[1][2] = sin(beta);
    result[2][1] = -sin(beta);
    result[2][2] = cos(beta);
    result[3][3] = 1;
}

void NumericalFlux::A_minus(double result[4][4], double w0, double w1, double w3, double w4)
{
    double _R[4][4];
    double _D_minus[4][4];
    double _R_inv[4][4];
    double _A_minus[4][4];
    double _tmp[4][4];
    for (int component=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _R[i][j] = matrix_R(i, j, w0, w1, w3, w4);
    for (int component=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _D_minus[i][j] = matrix_D_minus(i, j, w0, w1, w3, w4);
    for (int component=0; i < 4; i++)
        for (int j=0; j < 4; j++)
            _R_inv[i][j] = matrix_R_inv(i, j, w0, w1, w3, w4);
    dot(_tmp, _D_minus, _R_inv);
    dot(result, _R, _tmp);
}

void NumericalFlux::riemann_solver(double result[4], double w_L[4], double w_R[4])
{
    //printf("w_l: %f %f %f %f\n", w_L[0], w_L[1], w_L[2], w_L[3]);
    //printf("w_r: %f %f %f %f\n", w_R[0], w_R[1], w_R[2], w_R[3]);
    double _tmp1[4][4];
    double _tmp2[4][4];
    double _tmp3[4];
    double _tmp4[4];
    A_minus(_tmp1, w_R[0], w_R[1], w_R[2], w_R[3]);
    A_minus(_tmp2, w_L[0], w_L[1], w_L[2], w_L[3]);
    dot_vector(_tmp3, _tmp1, w_r);
    dot_vector(_tmp4, _tmp2, w_l);
    for (int component=0; i < 4; i++) {
        double _1 = f_x(i, w_L[0], w_L[1], w_L[2], w_L[3]);
        double _2 = _tmp3[i];
        double _3 = _tmp4[i];
        result[i] = _1 + _2 - _3;
    }
}

// calculates the iget_nvert()ed flux, for testing purposes
// it should return the same thing as riemann_solver(), only with minus sign
void NumericalFlux::riemann_solver_iget_nvert()(double result[4], double w_L[4], double w_R[4])
{
    double m[4][4];
    double _w_L[4];
    double _w_R[4];
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
// state as a 4-component vector.
void NumericalFlux::numerical_flux(double result[4], double w_L[4], double w_R[4],
        double nx, double ny)
{
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
}

// The same as numerical_flux, but only returns the i-th component:
double NumericalFlux::numerical_flux_i(int component, double w_L[4], double w_R[4],
        double nx, double ny)
{
    double result[4];
    numerical_flux(result, w_l, w_r, nx, ny);
    return result[i];
}
*/