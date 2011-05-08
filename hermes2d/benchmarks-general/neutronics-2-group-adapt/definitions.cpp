#include "weakform_library/weakforms_neutronics.h"
#include "function/filter.h"

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
                   const std::vector<DefaultFunction*>& f_src)
      : DefaultWeakFormFixedSource(matprop, f_src)
    {
      for (unsigned int g = 0; g < matprop.get_G(); g++)
      {
        add_matrix_form_surf(new GammaBoundaryCondition::Jacobian(g, bdy_gamma, matprop));
        add_vector_form_surf(new GammaBoundaryCondition::Residual(g, bdy_gamma, matprop));
      }
    }
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
    
    std::string get_material(double x, double y) const
    {
      if (x >= a && x < c && y >= a && y < c) return regions[0];
      if (x > c && x <= b && y >= a && y < c) return regions[1];
      if (x > c && x <= b && y > c && y <= b) return regions[2];
      if (x >= a && x < c && y > c && y <= b) return regions[3];
      return std::string();
    }
};

class CustomRightHandSide_g1 : public DefaultFunction, public CustomPiecewiseFunction
{
public:
  CustomRightHandSide_g1(double a, double b, const MaterialPropertyMaps& matprop, const std::string regions[])
    : DefaultFunction(), CustomPiecewiseFunction(a, b, matprop, regions)
  {};

  virtual double value(double x, double y) const {
    std::string mat = get_material(x,y);
        
    rank1 nSf = WeakFormsNeutronics::Multigroup::MaterialProperties::Common::NDArrayMapOp::multiply<rank1>(
      matprop.get_nu(mat), matprop.get_Sigma_f(mat)
    );
    
    rank1 D = matprop.get_D(mat);
    rank1 Sr = matprop.get_Sigma_r(mat);

    double L =  1./40.*exp(-4*sqr(x)) * ( 
                  20*( 1+4*(8*sqr(x)-1)*y*(y-2) ) * D[0] + 1./8.*y*(y-2) * (
                    80*nSf[0] + (
                      10-2*cos(8*M_PI*x) + cos(8*M_PI*(x-y)) - 2*cos(8*M_PI*y) + cos(8*M_PI*(x+y))
                    )*nSf[1] - 80*Sr[0] 
                  )
                );
    return L;
  }

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

  virtual double value(double x, double y) const {
    std::string mat = get_material(x,y);
    
    rank1 D = matprop.get_D(mat);
    rank1 Sr = matprop.get_Sigma_r(mat);
    rank2 Ss = matprop.get_Sigma_s(mat);

    double yym2 = (y-2)*y;
    double p2 = sqr(M_PI), x2 = sqr(x), px = M_PI*x, py = M_PI*y;
    double s8px = sin(8*px);
    double s8py = sin(8*py);
    double c8px = cos(8*px);
    double c8py = cos(8*py);
    double s4px2 = sqr(sin(4*px));
    double s4py2 = sqr(sin(4*py));

    double L =  1./40.*exp(-4*x2) * (
                  10*yym2*Ss[1][0] - yym2*Sr[1]*( 1+s4px2*s4py2 ) + 0.5*D[1]*(
                    5 + 20*(8*x2-1)*yym2 + (
                      -1 + 4*(1+8*p2-8*x2)*yym2
                    )*c8py + c8px*( 
                      -1 + 4*(1+8*p2-8*x2)*yym2 + ( 1 - 4*(1+16*p2-8*x2)*yym2 )*c8py
                    ) + 32*M_PI*(
                      -4*x*yym2*s8px*s4py2 + (y-1)*s4px2*s8py
                    )
                  )
                );  
    return L; 
  }

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

  double value(double x, double y) const 
  {
    return exp(-4*sqr(x))*(y/2.-sqr(y/2.));
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    double em4x2 = exp(-4*sqr(x));
    dx =  2.0*em4x2*x*y*(y-2);
    dy = -0.5*em4x2*(y-1);
  };

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(Ord::get_max_order());
  }
};

class CustomExactSolution_g2 : public ExactSolutionScalar
{
public:
  CustomExactSolution_g2(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const 
  {
    double px = M_PI*x, py = M_PI*y;
    double x2 = sqr(x);
    double em4x2 = exp(-4*x2);
    double s4px2 = sqr(sin(4*px));
    double s4py2 = sqr(sin(4*py));
    return em4x2*(y/2.-sqr(y/2.)) * (1 + s4px2 * s4py2) / 10.0;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    double yym2 = (y-2)*y;
    double x2 = sqr(x), px = M_PI*x, py = M_PI*y;
    double em4x2 = exp(-4*x2);
    double s8px = sin(8*px);
    double s8py = sin(8*py);
    double s4px2 = sqr(sin(4*px));
    double s4py2 = sqr(sin(4*py));
    
    dx =  1./10.*em4x2*yym2 * ( 2*x + (2*x*s4px2 - M_PI*s8px) * s4py2 );
    dy =  1./20.*em4x2 * ( 1-y-s4px2*( (y-1)*s4py2 + 2*M_PI*yym2*s8py )  );
  };

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(Ord::get_max_order());
  }
};

//////  Physical parameters.  /////////////////////////////////////////////////////////////////
// Two-group material properties for the 4 macro regions.

const std::string regions[4] = {
  "region1", "region2", "region3", "region4"
};

/*
const MaterialPropertyMap1 D = material_property_map<rank1>
(
  regions[0], grow(1.12)(0.6)
)(
  regions[1], grow(1.2)(0.5)
)(
  regions[2], grow(1.35)(0.8)
)(
  regions[3], grow(1.3)(0.9)
);
*/
// FIXME: The benchmark and its exact solution must be fixed for discontinuous D

const MaterialPropertyMap1 D = material_property_map<rank1>
(
  regions[0], grow(1)(0.5)
)(
  regions[1], grow(1)(0.5)
)(
  regions[2], grow(1)(0.5)
)(
  regions[3], grow(1)(0.5)
);

const MaterialPropertyMap1 Sr = material_property_map<rank1>
(
  regions[0], grow(0.011)(0.13)
)(
  regions[1], grow(0.09)(0.15)
)(
  regions[2], grow(0.035)(0.25)
)(
  regions[3], grow(0.04)(0.35)
);

const MaterialPropertyMap1 nSf = material_property_map<rank1>
(
  regions[0], grow(0.0025)(0.15)
)(
  regions[1], grow(0.00)(0.00)
)(
  regions[2], grow(0.0011)(0.1)
)(
  regions[3], grow(0.004)(0.25)
);

const double nu = 2.43;

const double chi_data[2] = {1, 0};
const rank1 chi(chi_data, chi_data+2);

const MaterialPropertyMap2 Ss = material_property_map<rank2>
(
  regions[0],
  gmat
  (
    grow(0.0)(0.0)
  )(
    grow(0.05)(0.0)
  )
)(
  regions[1],
  gmat
  (
    grow(0.0)(0.0)
  )(
    grow(0.08)(0.0)
  )
)(
  regions[2],
  gmat
  (
    grow(0.0)(0.0)
  )(
    grow(0.025)(0.0)
  )
)(
  regions[3],
  gmat
  (
    grow(0.0)(0.0)
  )(
    grow(0.014)(0.0)
  )
);

const bool2 scattering_mg_structure = bool2
(
  bool_mat
  (
    bool_row(false)(false)
  )(
    bool_row(true)(false)
  )
);
const bool1 fission_mg_structure = bool1
(
  bool_row(true)(false)
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Other utility functions:

// Construct a string representing the globally set adaptivity options.
void make_str_from_adapt_opts(std::stringstream& str)
{    
  switch (CAND_LIST) {
    case H2D_H_ANISO:
    case H2D_H_ISO:
      str << "h" << P_INIT[0] << "_" << P_INIT[1];
      break;
    case H2D_P_ANISO:
    case H2D_P_ISO:
      str << "p" << INIT_REF_NUM[0] << "_" << INIT_REF_NUM[1];
      break;
    default:
      str << "hp";
      break;
  }
  switch (CAND_LIST) {
    case H2D_H_ANISO:
    case H2D_P_ANISO:
    case H2D_HP_ANISO:
      str << "_aniso";
      break;
    case H2D_H_ISO:
    case H2D_P_ISO:
    case H2D_HP_ISO:
      str << "_iso";
      break;
    case H2D_HP_ANISO_H:
      str << "_anisoh";
      break;
    case H2D_HP_ANISO_P:
      str << "_anisop";
      break;
  }
  
  str << (MULTIMESH ? "_multi" : "_single");
}


