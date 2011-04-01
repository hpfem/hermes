#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"

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

class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(Hermes::vector<Space *> spaces, 
                        Hermes::vector<Solution *> solutions);

  /// Destructor.
   ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids(double threshold);

  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  double calculate_jumps(Element* e, int edge_i);

  /// Calculates the norm of the solution on the central element.
  double calculate_norm(Element* e, int edge_i);

protected:
  /// Members.
  Hermes::vector<Space *> spaces;
  Hermes::vector<Solution *> solutions;
  Mesh* mesh;
  std::set<int> discontinuous_element_ids;
};

class FluxLimiter
{
public:
  /// Constructor.
  FluxLimiter(scalar* solution_vector, Hermes::vector<Space *> spaces, Hermes::vector<Solution *> solutions);

  /// Destructor.
   ~FluxLimiter();

  /// Do the limiting.
  void limit_according_to_detector(std::set<int>& discontinuous_elements);

protected:
  /// Members.
  scalar* solution_vector;
  Hermes::vector<Space *> spaces;
  Hermes::vector<Solution *> solutions;
};

// Filters.
class MachNumberFilter : public SimpleFilter
{
public: 
  MachNumberFilter(Hermes::vector<MeshFunction*> solutions, double kappa) : SimpleFilter(solutions), kappa(kappa) {};
  ~MachNumberFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);

  double kappa;
};

class PressureFilter : public SimpleFilter
{
public: 
  PressureFilter(Hermes::vector<MeshFunction*> solutions, double kappa) : SimpleFilter(solutions), kappa(kappa) {};
  ~PressureFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);

  double kappa;
};

class EntropyFilter : public SimpleFilter
{
public: 
  EntropyFilter(Hermes::vector<MeshFunction*> solutions, double kappa, double rho_ext, double p_ext) : SimpleFilter(solutions), kappa(kappa), rho_ext(rho_ext), p_ext(p_ext) {};
  ~EntropyFilter() {};
protected:
  virtual void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result);

  double kappa, rho_ext, p_ext;
};

#endif