#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace RefinementSelectors;
using namespace Views;

#define LINEAR_NONLINEAR_SWITCH

class CustomNonlinearity : public Hermes1DFunction<double>
{
public:
  CustomNonlinearity(double alpha);

  virtual double value(double u) const;

  virtual Ord value(Ord u) const;

  virtual double derivative(double u) const;

  virtual Ord derivative(Ord u) const;

  protected:
    double alpha;
};

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  CustomEssentialBCNonConst(std::string marker);

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
  ~CustomInitialCondition();

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const;
};

class CustomWeakFormSteadyState : public WeakForm<double>
{
public:
  CustomWeakFormSteadyState(Hermes1DFunction<double>* thermal_conductivity, Hermes2DFunction<double>* heat_source);
};
  
class CustomWeakFormTimeDependent : public CustomWeakFormSteadyState
{
public:
  CustomWeakFormTimeDependent(Hermes1DFunction<double>* thermal_conductivity, Hermes2DFunction<double>* heat_source, MeshFunctionSharedPtr<double> prev_sln);
  
  class CustomMatrixFormVol : public DefaultMatrixFormVol<double>
  {
    CustomMatrixFormVol(int i, int j);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual MatrixFormVol<double>* clone() const;

    friend class CustomWeakFormTimeDependent;
  };
  
  class CustomVectorFormVol : public DefaultResidualVol<double>
  {
    CustomVectorFormVol(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual VectorFormVol<double>* clone() const;

    friend class CustomWeakFormTimeDependent;
  };
  
  class CustomResidualFormVol : public DefaultResidualVol<double>
  {
    CustomResidualFormVol(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual VectorFormVol<double>* clone() const;

    friend class CustomWeakFormTimeDependent;
  };
};

void load_mesh(MeshSharedPtr& mesh, const char* filename, int num_initial_refinements);

void display(MeshFunctionSharedPtr<double>& sln_to_display, SpaceSharedPtr<double>& space_to_display);

double* get_initial_Newton_guess(int as, WeakForm<double>* wf, SpaceSharedPtr<double> space, SpaceSharedPtr<double> ref_space, MeshFunctionSharedPtr<double> sln);

#define HERMES_ANY_MARKER "-1234"