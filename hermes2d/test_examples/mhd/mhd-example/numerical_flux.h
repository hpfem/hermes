#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H
#include "util.h"

class NumericalFlux
{
public:
  NumericalFlux(double kappa);

  /// Calculates a specified component of the flux on an interior edge.
  /// \param[in] component
  /// \param[in] w_L The Left ("central") state
  /// \param[in] w_L The Right ("neighbor") state
  /// \param[in] nx, ny Components of the unit outer normal vector
  /// Returns the result.
  virtual double numerical_flux(int component, double w_L[4], double w_R[4], double nx, double ny) = 0;

  /// Calculates a specified component of the flux on the inlet.
  /// \param[in] component
  /// \param[in] w_L The Left ("central") state
  /// \param[in] w_B The boundary state
  /// \param[in] nx, ny Components of the unit outer normal vector
  /// Returns the result.
  virtual double numerical_flux_inlet(int component, double w_L[4], double w_B[4], double nx, double ny) = 0;

  /// Calculates a specified component of the flux on the outlet.
  /// \param[in] component
  /// \param[in] w_L The Left ("central") state
  /// \param[in] w_B The boundary state
  /// \param[in] nx, ny Components of the unit outer normal vector
  /// Returns the result.
  virtual double numerical_flux_outlet(int component, double w_L[4], double w_B[4], double nx, double ny) = 0;

  /// Rotates the state_vector into the local coordinate system.
  /// Used especially in solid wall boundary condition - the state on the boundary is the inner state with velocity reflected off the wall.
  void Q(double result[4], double state_vector[4], double nx, double ny);

  /// Rotates the state_vector back from the local coordinate system.
  /// Used especially in solid wall boundary condition - the state on the boundary is the inner state with velocity reflected off the wall.
  void Q_inv(double result[4], double state_vector[4], double nx, double ny);

  // Poisson adiabatic constant = c_p/c_v = 1 + R/c_v.
  double kappa;
};

class UpWindNumericalFlux : public NumericalFlux
{
public:
  UpWindNumericalFlux(double kappa);

  virtual double numerical_flux(int component, double w_L[4], double w_R[4], double nx, double ny);

  virtual double numerical_flux_inlet(int component, double w_L[4], double w_B[4], double nx, double ny);
  
  virtual double numerical_flux_outlet(int component, double w_L[4], double w_B[4], double nx, double ny);
};

#endif
