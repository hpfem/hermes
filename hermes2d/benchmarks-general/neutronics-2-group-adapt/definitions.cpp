#include "definitions.h"
   
CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop, std::string bdy_gamma,
                                const std::vector<DefaultFunction*>& f_src )
  : DefaultWeakFormFixedSource(matprop, f_src)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_matrix_form_surf(new GammaBoundaryCondition::Jacobian(g, bdy_gamma, matprop));
    add_vector_form_surf(new GammaBoundaryCondition::Residual(g, bdy_gamma, matprop));
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Right hand sides:
    
std::string CustomPiecewiseFunction::get_material(double x, double y) const
{
  if (x >= a && x < c && y >= a && y < c) return regions[0];
  if (x > c && x <= b && y >= a && y < c) return regions[1];
  if (x > c && x <= b && y > c && y <= b) return regions[2];
  if (x >= a && x < c && y > c && y <= b) return regions[3];
  return std::string();
}

double CustomRightHandSide_g1::value(double x, double y) const 
{
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


double CustomRightHandSide_g2::value(double x, double y) const 
{
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution:

inline double CustomExactSolution_g1::value(double x, double y) const 
{
  return exp(-4*sqr(x))*(y/2.-sqr(y/2.));
}

void CustomExactSolution_g1::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  double em4x2 = exp(-4*sqr(x));
  dx =  2.0*em4x2*x*y*(y-2);
  dy = -0.5*em4x2*(y-1);
}

double CustomExactSolution_g2::value(double x, double y) const 
{
  double px = M_PI*x, py = M_PI*y;
  double x2 = sqr(x);
  double em4x2 = exp(-4*x2);
  double s4px2 = sqr(sin(4*px));
  double s4py2 = sqr(sin(4*py));
  return em4x2*(y/2.-sqr(y/2.)) * (1 + s4px2 * s4py2) / 10.0;
}

void CustomExactSolution_g2::derivatives (double x, double y, scalar& dx, scalar& dy) const 
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
}


