#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"
#include "discrete_problem/dg/discrete_problem_dg_assembler.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

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

class CFLCalculation
{
public:
  CFLCalculation(double CFL_number, double kappa);

  // If the time step is necessary to decrease / possible to increase, the value time_step will be rewritten.
  bool calculate(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step) const;
  bool calculate(double* sln_vector, Hermes::vector<SpaceSharedPtr<double> > spaces, double & time_step) const;
  void calculate_semi_implicit(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step) const;

  void set_number(double new_CFL_number);

protected:
  double CFL_number;
  double kappa;
};

class ADEStabilityCalculation
{
public:
  ADEStabilityCalculation(double AdvectionRelativeConstant, double DiffusionRelativeConstant, double epsilon);

  // If the time step is necessary to decrease / possible to increase, the value time_step will be rewritten.
  void calculate(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, MeshSharedPtr mesh, double & time_step);

protected:
  double AdvectionRelativeConstant;
  double DiffusionRelativeConstant;
  double epsilon;
};

class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, 
    Hermes::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
  ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  virtual std::set<int>& get_discontinuous_element_ids() = 0;
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRho() { return this->oscillatory_element_idsRho; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoVX() { return this->oscillatory_element_idsRhoVX; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoVY() { return this->oscillatory_element_idsRhoVY; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoE() { return this->oscillatory_element_idsRhoE; }

protected:
  /// Members.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions;
  Hermes::vector<Solution<double>*> solutionsInternal;
  std::set<int> discontinuous_element_ids;
  std::set<std::pair<int, double> > oscillatory_element_idsRho;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVX;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVY;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoE;
  MeshSharedPtr mesh;
};

class KrivodonovaDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KrivodonovaDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, 
    Hermes::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
  ~KrivodonovaDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();
  std::set<int>& get_discontinuous_element_ids(double threshold);

protected:
  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  void calculate_jumps(Element* e, int edge_i, double result[4]);

  /// Calculates h.
  double calculate_h(Element* e, int polynomial_order);

  /// Calculates the norm of the solution on the central element.
  void calculate_norms(Element* e, int edge_i, double result[4]);
};

class KuzminDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KuzminDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, 
    Hermes::vector<MeshFunctionSharedPtr<double> > solutions, bool limit_all_orders_independently = false);

  /// Destructor.
  ~KuzminDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();

  /// Return a reference to the inner structures.
  std::set<int>& get_second_order_discontinuous_element_ids();

  /// Returns info about the method.
  bool get_limit_all_orders_independently();
protected:
  /// Center.
  void find_centroid_values(Hermes::Hermes2D::Element* e, double u_c[4], double x_ref, double y_ref);
  void find_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dx_c[4], double u_dy_c[4], double x_ref, double y_ref);
  void find_second_centroid_derivatives(Hermes::Hermes2D::Element* e, double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4]);

  /// Vertices.
  void find_vertex_values(Hermes::Hermes2D::Element* e, double vertex_values[4][4]);
  void find_vertex_derivatives(Hermes::Hermes2D::Element* e, double vertex_derivatives[4][4][2]);

  /// Logic - 1st order.
  void find_u_i_min_max_first_order(Hermes::Hermes2D::Element* e, double u_i_min[4][4], double u_i_max[4][4]);
  void find_alpha_i_first_order(double u_i_min[4][4], double u_i_max[4][4], double u_c[4], double u_i[4][4], double alpha_i[4]);
  void find_alpha_i_first_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4], double u_c[4], double u_dx_c[4], double u_dy_c[4], double alpha_i_real[4]);

  /// Logic - 2nd order.
  void find_u_i_min_max_second_order(Hermes::Hermes2D::Element* e, double u_d_i_min[4][4][2], double u_d_i_max[4][4][2]);
  void find_alpha_i_second_order(double u_d_i_min[4][4][2], double u_d_i_max[4][4][2], double*** u_dx_c, double*** u_dy_c, double u_d_i[4][4][2], double alpha_i[4]);
  void find_alpha_i_second_order_real(Hermes::Hermes2D::Element* e, double u_i[4][4][2], double u_dx_c[4], double u_dy_c[4], double u_dxx_c[4], double u_dxy_c[4], double u_dyy_c[4], double alpha_i_real[4]);

private:
  /// For limiting of second order terms.
  std::set<int> second_order_discontinuous_element_ids;
  bool limit_all_orders_independently;
};

class FluxLimiter
{
public:
  /// Enumeration of types.
  /// Used to pick the proper DiscontinuityDetector.
  enum LimitingType
  {
    Krivodonova,
    Kuzmin
  };
  /// Constructor.
  FluxLimiter(LimitingType type, double* solution_vector, Hermes::vector<SpaceSharedPtr<double> > spaces, bool Kuzmin_limit_all_orders_independently = false);
  FluxLimiter(FluxLimiter::LimitingType type, Hermes::vector<MeshFunctionSharedPtr<double> > solutions, Hermes::vector<SpaceSharedPtr<double>  > spaces, bool Kuzmin_limit_all_orders_independently = false);

  /// Destructor.
  ~FluxLimiter();

  /// Do the limiting.
  /// With the possibility to also limit the spaces from which the spaces in the constructors are refined.
  virtual int limit_according_to_detector(Hermes::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit = Hermes::vector<SpaceSharedPtr<double> >());

  /// For Kuzmin's detector.
  virtual void limit_second_orders_according_to_detector(Hermes::vector<SpaceSharedPtr<double> > coarse_spaces_to_limit = Hermes::vector<SpaceSharedPtr<double> >());

  void get_limited_solutions(Hermes::vector<MeshFunctionSharedPtr<double> > solutions_to_limit);
  bool limitOscillations;
protected:
  /// Members.
  double* solution_vector;
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  DiscontinuityDetector* detector;
  Hermes::vector<MeshFunctionSharedPtr<double> > limited_solutions;
};

// Filters.
class MachNumberFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
  ~MachNumberFilter() 
  {
  };

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }

    MachNumberFilter* filter = new MachNumberFilter(slns, this->kappa);

    return filter;
  }

protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
};

class PressureFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
  ~PressureFilter() 
  {
  };

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    PressureFilter* filter = new PressureFilter(slns, this->kappa);

    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
};

class VelocityFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public:
  // Vector of solutions: 0-th position - density, 1-st position - velocity component.
  VelocityFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : SimpleFilter<double>(solutions) {};
  ~VelocityFilter() 
  {
  };

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }

    VelocityFilter* filter = new VelocityFilter(slns);
    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
};

class EntropyFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  EntropyFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa, double rho_ext, double p_ext) : SimpleFilter<double>(solutions), kappa(kappa), rho_ext(rho_ext), p_ext(p_ext) {};
  ~EntropyFilter() 
  {
  };
  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
      slns.push_back(this->sln[i]->clone());
    EntropyFilter* filter = new EntropyFilter(slns, this->kappa, rho_ext, p_ext);

    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa, rho_ext, p_ext;
};

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

enum EulerLimiterType
{
  VertexBased = 0,
  VertexBasedWithLimitingNonConservative = 1,
  CoarseningJumpIndicatorDensity = 2,
  CoarseningJumpIndicatorDensityToAll = 3,
  CoarseningJumpIndicatorAllToThemselves = 4,
  CoarseningJumpIndicatorAllToAll = 5,
  VertexBasedPCoarsener = 6
};

template<typename LimiterType>
void limitVelocityAndEnergy(Hermes::vector<SpaceSharedPtr<double> > spaces, LimiterType* limiter, Hermes::vector<MeshFunctionSharedPtr<double> > slns);

class FeistauerJumpDetector : public PostProcessing::Limiter<double>
{
public:
  FeistauerJumpDetector(SpaceSharedPtr<double> space, double* solution_vector);
  FeistauerJumpDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, double* solution_vector);
  void init();
  ~FeistauerJumpDetector();

  void set_type(EulerLimiterType indicatorType);
  EulerLimiterType get_type();

  /// Alpha in the indicator, it is the exponent (h^{alpha}) in the denominator.
  /// Should be between 1 - 5.
  static double alpha;

  /// The constant the jump indicator g is compared to and when larger than this, the limiting is carried out.
  static double thresholdConstant;

private:
  void process();
  void get_jump_indicators(Element* e, double* values);
  void assemble_one_neighbor(NeighborSearch<double>& ns, int edge, unsigned int neighbor_i, double* values);
  bool conditionally_coarsen(double max_value, double* values, Element* e);
  
  EulerLimiterType indicatorType;
};

class DensityErrorCalculator
{
public:
  DensityErrorCalculator();
  ~DensityErrorCalculator();
  void process(MeshFunctionSharedPtr<double> density, SpaceSharedPtr<double> density_space);

  double* element_errors;
  double4* edge_errors;
  
private:
  void init();
  void get_element_error(Element* e, int order);
  void get_edges_error(Element* e, int order);
  void assemble_one_neighbor(NeighborSearch<double>& ns, int edge, unsigned int neighbor_i, int order);

  // number of elements * (1 + number of edges) - for the volumetric and jump errors.
  int element_alloc_size;
  MeshFunctionSharedPtr<double> density;
};

PostProcessing::Limiter<double>* create_limiter(EulerLimiterType limiter_type, SpaceSharedPtr<double> space, double* solution_vector, int polynomial_degree = 1, bool verbose = false);
PostProcessing::Limiter<double>* create_limiter(EulerLimiterType limiter_type, Hermes::vector<SpaceSharedPtr<double> > spaces, double* solution_vector, int polynomial_degree = 1, bool verbose = false);

#endif
