#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion; 

class CustomWeakForm : public DefaultWeakFormFixedSource
{
  struct GammaBoundaryCondition
  {
    static const double gamma = 8;
    
    class Jacobian : public WeakForm::MatrixFormSurf
    {
      public:
        Jacobian(unsigned int g, std::string bdy_gamma, const MaterialPropertyMaps& matprop) 
          : WeakForm::MatrixFormSurf(g,g,bdy_gamma), 
            g(g), matprop(matprop)
        {};
        
        scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
        {
          std::string mat = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
          return gamma * matprop.get_D(mat)[g] * int_u_v<double, scalar>(n, wt, u, v);
        }

        Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
        {
          return int_u_v<Ord, Ord>(n, wt, u, v);
        }
                      
      private:
        unsigned int g;
        const MaterialPropertyMaps& matprop;
    };
    
    class Residual : public WeakForm::VectorFormSurf
    {
      public:
        Residual(unsigned int g, std::string bdy_gamma, const MaterialPropertyMaps& matprop) 
          : WeakForm::VectorFormSurf(g,bdy_gamma), 
            g(g), matprop(matprop) 
        {};
        
        scalar value(int n, double *wt, Func<scalar> *u_ext[],
                     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
        {
          std::string mat = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
          return gamma * matprop.get_D(mat)[g] * int_u_v<double, scalar>(n, wt, u_ext[g], v);
        }

        Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
        {
          return int_u_v<Ord, Ord>(n, wt, u_ext[g], v);
        }
                      
      private:
        unsigned int g;
        const MaterialPropertyMaps& matprop;
    };
  };
  
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop, std::string bdy_gamma,
                   const std::vector<DefaultFunction*>& f_src);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Right hand sides:

class CustomPiecewiseFunction
{
  private:
    double a, b, c;    
    
  protected:
    const MaterialPropertyMaps& matprop;
    const std::string* regions;

  public:  
    CustomPiecewiseFunction(double a, double b, const MaterialPropertyMaps& matprop, const std::string regions[])
      : a(a), b(b), c((a+b)/2.), matprop(matprop), regions(regions)
    {};
    
    std::string get_material(double x, double y) const;
};

class CustomRightHandSide_g1 : public DefaultFunction, public CustomPiecewiseFunction
{
public:
  CustomRightHandSide_g1(double a, double b, const MaterialPropertyMaps& matprop, const std::string regions[])
    : DefaultFunction(), CustomPiecewiseFunction(a, b, matprop, regions)
  {};

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }
};

class CustomRightHandSide_g2 : public DefaultFunction, public CustomPiecewiseFunction
{
public:
  CustomRightHandSide_g2(double a, double b, const MaterialPropertyMaps& matprop, const std::string regions[])
    : DefaultFunction(), CustomPiecewiseFunction(a, b, matprop, regions)
  {};

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution:

class CustomExactSolution_g1 : public ExactSolutionScalar
{
public:
  CustomExactSolution_g1(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const;
  
  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(Ord::get_max_order());
  }
};

class CustomExactSolution_g2 : public ExactSolutionScalar
{
public:
  CustomExactSolution_g2(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const;
  
  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(Ord::get_max_order());
  }
};
