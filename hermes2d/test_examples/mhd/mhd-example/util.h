#ifndef _UTIL_H
#define _UTIL_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// Class calculating various quantities
class QuantityCalculator
{
public:
  // Calculates energy from other quantities.
  static double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa);
 
  // Calculates pressure from other quantities.
  static double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
 
  // Calculates speed of sound.
  static double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
};

///TODO Tato trida by se mela rozsirit o zbyle slozky fluxu a take by se stavajici fluxy meli upravit pokud je treba
class EulerFluxes
{
public:
  EulerFluxes(double kappa) : kappa(kappa) {}

  
  double A_1_0_0(double rho, double rho_v_x, double rho_v_y, double energy) {

    return double(0.0);
  }

  
  double A_1_0_1(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(1.0);
  }

  
  double A_1_0_2(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(0.0);
  }

  
  double A_1_0_3(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(0.0);
  }

  
  double A_2_0_0(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(0.0);
  }

  
  double A_2_0_1(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(0.0);
  }

  
  double A_2_0_2(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(1.0);
  }

  
  double A_2_0_3(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(0.0);
  }

  
  double A_1_1_0(double rho, double rho_v_x, double rho_v_y, double energy) {
    return double(- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (kappa - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  double A_1_1_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((3. - kappa) * (rho_v_x / rho));
  }

  
  double A_1_1_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((1.0 - kappa) * (rho_v_y / rho));
  }

  
  double A_1_1_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(kappa - 1.);
  }

  
  double A_2_1_0(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(- rho_v_x * rho_v_y / (rho * rho));
  }

  
  double A_2_1_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(rho_v_y / rho);
  }

  
  double A_2_1_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(rho_v_x / rho);
  }

  
  double A_2_1_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(0);
  }

  
  double A_1_2_0(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(- rho_v_x * rho_v_y / (rho * rho));
  }

  
  double A_1_2_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(rho_v_y / rho);
  }

  
  double A_1_2_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(rho_v_x / rho);
  }

  
  double A_1_2_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(0);
  }

  
  double A_2_2_0(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (kappa - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  double A_2_2_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((1.0 - kappa) * (rho_v_x / rho));
  }

  
  double A_2_2_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((3.0 - kappa) * (rho_v_y / rho));
  }

  
  double A_2_2_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(kappa - 1.);
  }

  
  double A_1_3_0(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((rho_v_x / rho) * (((kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (kappa * energy / rho)));
  }

  
  double A_1_3_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((kappa * energy / rho) - (kappa - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (kappa - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  
  double A_1_3_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  double A_1_3_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(kappa * (rho_v_x / rho));
  }

  
  double A_2_3_0(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));
  }

  
  double A_2_3_1(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  double A_2_3_2(double rho, double rho_v_x, double rho_v_y, double energy){
    return double((energy / rho) + (1 / rho) * (kappa - 1.0) * ( energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - kappa) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  
  double A_2_3_3(double rho, double rho_v_x, double rho_v_y, double energy){
    return double(kappa * rho_v_y / rho);
  }
  protected:
    double kappa;
};

#endif
