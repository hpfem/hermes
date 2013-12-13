#include "definitions.h"
#include "forms.h"

char* SolvedExampleString[5] = { "1D", "CircularConvection", "MovingPeak", "AdvectedCube", "SolidBodyRotation" };

double upwind_flux(double u_cent, double u_neib, double a_dot_n)
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord upwind_flux(Ord u_cent, Ord u_neib, double a_dot_n)
{
  return a_dot_n * (u_cent + u_neib);
}

static double advection_term_cube(double x, double y, double vx, double vy)
{
  return vx + vy;
}

static double advection_term_solid_body_rotation(double x, double y, double vx, double vy)
{
  return vx * (0.5 - y) + vy * (x - 0.5);
}

static double advection_term_circular_convection(double x, double y, double vx, double vy)
{
  return (vx * y) + (vy * (1 - x));
}

static double advection_term_moving_peak(double x, double y, double vx, double vy)
{
  return -(vx * y) + (vy * x);
}

static double advection_term_benchmark(double x, double y, double vx, double vy)
{
  return vx + vy;
}
scalar_product_with_advection_direction advection_term;
bool only_x_der;

ErrorWeakForm::ErrorWeakForm(SolvedExample solvedExample)
{
  this->add_vector_form(new ErrorFormVol());
}

static void initialization(SolvedExample solvedExample)
{
  only_x_der = false;
  switch(solvedExample)
  {
  case AdvectedCube:
    advection_term = advection_term_cube;
    break;
  case SolidBodyRotation:
    advection_term = advection_term_solid_body_rotation;
    break;
  case CircularConvection:
    advection_term = advection_term_circular_convection;
    break;
  case MovingPeak:
    advection_term = advection_term_moving_peak;
    break;
  case Benchmark:
    advection_term = advection_term_benchmark;
    only_x_der = true;
    break;
  }
}

MassWeakForm::MassWeakForm() : WeakForm<double>(1) 
{
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
}

SmoothingWeakForm::SmoothingWeakForm(SolvedExample solvedExample, bool local, int explicitSchemeStep, bool add_inlet, std::string inlet , double diffusivity, double s, double sigma, bool add_rhs) : WeakForm<double>(1) 
{
  initialization(solvedExample);

  // Matrix
  // M
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  // A_tilde  
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, local));
  // A_tilde_surf
  this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));

  if(std::abs(diffusivity) > 1e-7)
  {
    add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
    add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, local, diffusivity, s, sigma));
    if(add_inlet)
      this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet));
  }

  // RHS
  // A
  add_vector_form(new CustomVectorFormVolConvection(0, 0));
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, true, true));

  if(std::abs(diffusivity) > 1e-7)
  {
    add_vector_form(new CustomVectorFormVolDiffusion(0, 0, diffusivity));
    add_vector_form_DG(new CustomVectorFormInterfaceDiffusion(0, 0, diffusivity, s, sigma));
  }

  // A_surf
  add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, false, true));

  if(std::abs(diffusivity) > 1e-7)
  {
    if(add_inlet)
      this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet, true));
  }
}

SmoothingWeakFormResidual::SmoothingWeakFormResidual(SolvedExample solvedExample, int explicitSchemeStep, bool add_inlet, std::string inlet , double diffusivity, double s, double sigma, bool add_rhs) : WeakForm<double>(1)
{
  initialization(solvedExample);

  // b
  if(add_inlet)
  {
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 1, true, false));
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 1, diffusivity, s, sigma, inlet, false, 1.));
  }

  add_vector_form(new CustomVectorFormVolConvection(0, 0));
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, true, true));
  if(std::abs(diffusivity) > 1e-7)
  {
    add_vector_form(new CustomVectorFormVolDiffusion(0, 0, diffusivity));
    add_vector_form_DG(new CustomVectorFormInterfaceDiffusion(0, 0, diffusivity, s, sigma));
  }

  // A_surf
  add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, false, true));

  if(std::abs(diffusivity) > 1e-7)
  {
    if(add_inlet)
      this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet, true));
  }
}

ExactWeakForm::ExactWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution) : WeakForm<double>(1) 
{
  initialization(solvedExample);

  // A_tilde  
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, false));
  // A_tilde_surf
  this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));

  if(std::abs(diffusivity) > 1e-7)
  {
    add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
    add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, false, diffusivity, s, sigma));
    if(add_inlet)
      this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet));
  }

  // b
  if(add_inlet)
  {
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, true, false));
    this->vfsurf.back()->set_ext(exact_solution);

    if(std::abs(diffusivity) > 1e-7)
    {
      this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet, false, 1.));
      this->vfsurf.back()->set_ext(exact_solution);
    }
  }
}

TimeDepWeakForm::TimeDepWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution) :   ExactWeakForm(solvedExample, add_inlet, inlet, diffusivity, s, sigma, exact_solution)
{
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
  add_vector_form(new CustomVectorFormVol(0, 0, 1.));
}

MultiscaleWeakForm::MultiscaleWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, bool local) : WeakForm<double>(1)
{
  initialization(solvedExample);

  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  // A_tilde  
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, local));
  // A_tilde_surf
  this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));

  if(std::abs(diffusivity) > 1e-7)
  {
    add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
    add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, local, diffusivity, s, sigma));
    if(add_inlet)
      this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet));
  }

  // b
  if(add_inlet)
  {
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, true, false));
    this->vfsurf.back()->set_ext(exact_solution);

    if(std::abs(diffusivity) > 1e-7)
    {
      this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet, false, 1.));
      this->vfsurf.back()->set_ext(exact_solution);
    }
  }
}

ExplicitWeakFormOffDiag::ExplicitWeakFormOffDiag(SolvedExample solvedExample, bool add_inlet, std::string inlet , double diffusivity, double s, double sigma) : WeakForm<double>(1) 
{
  initialization(solvedExample);

  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, true, true));

  if(std::abs(diffusivity) > 1e-7)
    add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, true, diffusivity, s, sigma, true));
}

ImplicitWeakForm::ImplicitWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet , double diffusivity, double s, double sigma) : WeakForm<double>(1)
{
  initialization(solvedExample);

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
  // Mass matrix - rhs
  add_vector_form(new CustomVectorFormVol(0, 0, 1.));

  // Convection.
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, false));
  // Numerical flux - inner-edge inlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, false, false));
  // Numerical flux - inner-edge inlet+outlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 1, true, true));
  // No convection - when test functions have zero gradient
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 2, true, false));
    // Numerical flux - boundary outlet - derivatives - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 1, false, true));
  }

  // Diffusion.
  // No Diffusion - when test functions have zero gradient
  add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
  add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, false, diffusivity, s, sigma));
  add_vector_form_DG(new CustomVectorFormInterfaceDiffusion(0, 1, diffusivity, s, sigma));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 1, diffusivity, s, sigma, inlet, true));
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 2, diffusivity, s, sigma, inlet, false, 1.));
  }
}

ExplicitWeakForm::ExplicitWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet , double diffusivity, double s, double sigma) : WeakForm<double>(1) 
{
  initialization(solvedExample);

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  // Convection.
  // Convective term - matrix
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  // Convective term - rhs
  add_vector_form(new CustomVectorFormVolConvection(0, 0));

  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, true));
  // Numerical flux - inner-edge inlet+outlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, true, true));
  // Numerical flux - inner-edge inlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 1, true, false));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 2, true, false));
    // Numerical flux - boundary outlet - mean values - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, false, true));
  }

  // Diffusion.
  // Diffusive term - matrix
  add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
  // Diffusive term - rhs - for mean values - zero
  add_vector_form(new CustomVectorFormVolDiffusion(0, 0, diffusivity));
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, false, diffusivity, s, sigma));
  add_vector_form_DG(new CustomVectorFormInterfaceDiffusion(0, 0, diffusivity, s, sigma));
  //add_vector_form_DG(new CustomVectorFormInterfaceDiffusionOffDiag(0, 1, diffusivity, s, sigma));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 0, diffusivity, s, sigma, inlet, true));
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 2, diffusivity, s, sigma, inlet, false, 1.));
  }

  // Mass matrix - rhs
  add_vector_form(new CustomVectorFormVol(0, 1, 1.));
}

InitialConditionAdvectedCube::InitialConditionAdvectedCube(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
{
}

void InitialConditionAdvectedCube::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;
  dy = 0.;
}

double InitialConditionAdvectedCube::value(double x, double y) const
{
  if(x < 0. && y < 0.0 && x > -1. && y > -1.)
    return 1.0;
  else
    return 0.0;
}

Ord InitialConditionAdvectedCube::ord(double x, double y) const 
{
  return Ord(1);
}

MeshFunction<double>* InitialConditionAdvectedCube::clone() const
{
  return new InitialConditionAdvectedCube(this->mesh);

}

void InitialConditionSolidBodyRotation::derivatives(double x, double y, double& dx, double& dy) const 
{

  double radius = 0.;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  {		
    dx = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
    dy = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
  }
  else
  {			
    //cone
    x_0 = 0.5;
    y_0 = 0.25;
    radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
    if((radius< 1.0)&&(x!=x_0)) 
    { 	
      dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
      dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
    }
    else
    {
      dx=0.; dy=0.;
    }	
  }
};

double InitialConditionSolidBodyRotation::value(double x, double y) const 
{

  double result = 0.0;
  double radius;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  { 
    result = (1.0+ std::cos(M_PI*radius))/4.0;
    return result;	
  }
  //slotted cylinder
  x_0 = 0.5;
  y_0 = 0.75;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius <= 1) 
  { 	
    if(fabs((x-x_0))>= 0.025) return 1.0;
    if(y>=0.85) return 1.0;
  }	
  //cone
  x_0 = 0.5;
  y_0 = 0.25;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius<= 1.0) 
  { 	
    result = 1.0-radius;
  }	
  return result;
};

Ord InitialConditionSolidBodyRotation::ord(double x, double y) const 
{
  return Ord(10);
};
MeshFunction<double>* InitialConditionSolidBodyRotation::clone() const
{
  return new InitialConditionSolidBodyRotation(this->mesh);
}


void ExactSolutionCircularConvection::derivatives(double x, double y, double& dx, double& dy) const 
{
  double radius = std::sqrt(std::pow(x - 1, 2.) + std::pow(y, 2.));
  double radius_dx = (2 * x - 2.) / radius;
  double radius_dy = 2 * y / radius;
  if(radius >= 0.5 && radius <= 0.8)
  {
    dx = -0.25 * std::sin(M_PI * ((radius) / 0.15)) * radius_dx;
    dx = -0.25 * std::sin(M_PI * ((radius) / 0.15)) * radius_dy;
  }
  else
    dx = dy = 0.;
};

double ExactSolutionCircularConvection::value(double x, double y) const 
{
  double radius = std::sqrt(std::pow(x - 1, 2.) + std::pow(y, 2.));
  if(radius >= 0.5 && radius <= 0.8)
    return 0.25 * (1. + std::cos(M_PI * ((radius - 0.65) / 0.15)));
  else
    return 0.;
};

Ord ExactSolutionCircularConvection::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionCircularConvection::clone() const
{
  return new ExactSolutionCircularConvection(this->mesh);
}

void ExactSolutionMovingPeak::set_current_time(double time)
{
  this->current_time = time;
  double x_0 = 0.;
  double y_0 = 0.5;
  this->x_hat = x_0 * std::cos(this->current_time) - y_0 * std::sin(this->current_time);
  this->y_hat = -x_0 * std::sin(this->current_time) + y_0 * std::cos(this->current_time);
}

void ExactSolutionMovingPeak::derivatives(double x, double y, double& dx, double& dy) const 
{
  double fraction_times_pi = 1. / (4. * this->diffusivity * this->current_time);
  double fn_value = this->value(x, y);
  dx = -fn_value * 2. * (x - this->x_hat) * fraction_times_pi;
  dy = -fn_value * 2. * (y - this->y_hat) * fraction_times_pi;
};

double ExactSolutionMovingPeak::value(double x, double y) const 
{
  double fraction_times_pi = 1. / (4. * this->diffusivity * this->current_time);
  double r_squared = ( (x - this->x_hat) * (x - this->x_hat) ) + ( (y - this->y_hat) * (y - this->y_hat) );
  return fraction_times_pi * std::pow(2.718281828, -r_squared * fraction_times_pi) / M_PI;
};

Ord ExactSolutionMovingPeak::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionMovingPeak::clone() const
{
  Solution<double>* newSln;

  if(this->get_type() == HERMES_SLN)
  {
    newSln = new Solution<double>;
    newSln->copy(this);
  }
  else
    newSln = new ExactSolutionMovingPeak(this->mesh, this->diffusivity, this->current_time);

  return newSln;
}

double ExactSolutionMovingPeak::get_current_time() const
{
  return this->current_time;
}



void InitialConditionBenchmark::derivatives(double x, double y, double& dx, double& dy) const 
{
  y = 0.;
  double fn_value = this->value(x, y);
  double denominator = 4. * (this->diffusivity * y + 0.001);
  dx = fn_value * (-2. * (x - 0.2) / denominator);

  double a = this->diffusivity;

  double fraction = -0.0141047 / std::pow(((a * y) + 0.001), 2.5);
  double exponenta = std::exp(-std::pow(x-y-0.2,2.) / (4*(a*y + 0.001)));

  double multiplier = (a*a*y) - (0.5*a*x*x) + (0.2*a*x) + (0.5*a*y*y) - (0.019*a) - (0.001*x) + (0.001*y) + 0.0002;
  dy = fraction * exponenta * multiplier;
};

double InitialConditionBenchmark::value(double x, double y) const 
{
  y = 0.;
  double fraction = 0.1 / (std::sqrt(4. * M_PI * (this->diffusivity * y + 0.001)));
  double denominator = 4. * (this->diffusivity * y + 0.001);
  double enumerator = -std::pow(x - y - 0.2, 2.);
  return fraction * std::exp(enumerator/denominator);
};

Ord InitialConditionBenchmark::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* InitialConditionBenchmark::clone() const
{
  Solution<double>* newSln;

  if(this->get_type() == HERMES_SLN)
  {
    newSln = new Solution<double>;
    newSln->copy(this);
  }
  else
    newSln = new InitialConditionBenchmark(this->mesh, this->diffusivity);

  return newSln;
}

void ExactSolutionBenchmark::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;
  double a = this->diffusivity;

  double fraction = -0.0141047 / std::pow(((a * y) + 0.001), 2.5);
  double exponenta = std::exp(-std::pow(x-y-0.2,2.) / (4*(a*y + 0.001)));

  double multiplier = (a*a*y) - (0.5*a*x*x) + (0.2*a*x) + (0.5*a*y*y) - (0.019*a) - (0.001*x) + (0.001*y) + 0.0002;
  dy = fraction * exponenta * multiplier;
};

double ExactSolutionBenchmark::value(double x, double y) const 
{
  double fraction = 0.1 / (std::sqrt(4. * M_PI * (this->diffusivity * y + 0.001)));
  double denominator = 4. * (this->diffusivity * y + 0.001);
  double enumerator = -std::pow(x - y - 0.2, 2.);
  return fraction * std::exp(enumerator/denominator);
};


Ord ExactSolutionBenchmark::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionBenchmark::clone() const
{
  return new ExactSolutionBenchmark(this->mesh, this->diffusivity);
}


void ExactSolutionBenchmark2::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;

  double a = this->diffusivity;

  double exponential = std::exp(-std::pow((5.*x)-(5.*y)-1.,2.)/(25.*((4.*a*y)+(l*l))));
  double a_1 = 2.*((5.*x)-(5.*y)-1.)*(2.*a*((5.*x)+(5.*y)-1.)+(5.*l*l));
  double a_2 = 35.*std::pow((4.*a*y)+(l*l), 2.) * std::sqrt(((4.*a*y) + (l*l)) / (l*l));
  double b_1 = -10.*a;
  double b_2 = 7.*l*l*std::pow(((4.*a*y)+(l*l))/(l*l), 3./2.);
  dy = exponential * ((a_1 / a_2) + (b_1 / b_2));
};

double ExactSolutionBenchmark2::value(double x, double y) const 
{
  double fraction = 5. / ( 7. * this->sigma(y));
  double denominator = this->l * this->sigma(y);
  double enumerator = x - y - (2./10.);
  return fraction * std::exp(-std::pow(enumerator/denominator, 2.));
};

double ExactSolutionBenchmark2::sigma(double y) const
{
  return std::sqrt(1. + ((4. * this->diffusivity * y) / (this->l * this->l)));
}

Ord ExactSolutionBenchmark2::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionBenchmark2::clone() const
{
  return new ExactSolutionBenchmark2(this->mesh, this->diffusivity);
}

double* merge_slns(double* solution_vector_coarse, SpaceSharedPtr<double> space_coarse, double* solution_vector_fine, SpaceSharedPtr<double> space_fine, SpaceSharedPtr<double> space_full, bool add)
{
  double* target = new double[space_full->get_num_dofs()];
  Element *e;
  for_all_active_elements(e, space_full->get_mesh())
  {
    AsmList<double> al_coarse, al_fine, al_target;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);
    space_full->get_element_assembly_list(e, &al_target);
    int shift = -1;
    if(space_fine->get_seq() == space_full->get_seq())
      shift = 0;

    for(int i = 0; i < al_coarse.cnt; i++)
      target[al_target.dof[i]] = solution_vector_coarse[al_coarse.dof[i]];

    if(add)
      for(int i = 0; i < al_coarse.cnt; i++)
        target[al_target.dof[i]] += solution_vector_fine[al_fine.dof[i]];

    for(int i = al_coarse.cnt; i < al_target.cnt; i++)
      target[al_target.dof[i]] = solution_vector_fine[al_fine.dof[i + shift]];
  }
  return target;
}

Hermes::Algebra::SimpleVector<double>* cut_off_linear_part(double* src_vector, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  SimpleVector<double>* vector = new SimpleVector<double>(space_coarse->get_num_dofs());
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    vector->set(al_coarse.dof[0],  src_vector[al_fine.dof[0]]);
  }
  return vector;
}

Hermes::Algebra::SimpleVector<double>* cut_off_quadratic_part(double* src_vector, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  SimpleVector<double>* vector = new SimpleVector<double>(space_coarse->get_num_dofs());
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    vector->set(al_coarse.dof[0],src_vector[al_fine.dof[0]]);
    vector->set(al_coarse.dof[1],src_vector[al_fine.dof[1]]);
    vector->set(al_coarse.dof[2],src_vector[al_fine.dof[2]]);
  }
  return vector;
}

Hermes::Algebra::SimpleVector<double>* cut_off_means(double* src_vector, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  SimpleVector<double>* vector = new SimpleVector<double>(space_coarse->get_num_dofs());
  vector->zero();
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    for(int i = 0; i < al_coarse.cnt; i++)
      vector->set(al_coarse.dof[i],  src_vector[al_fine.dof[i+1]]);
  }
  return vector;
}

void add_means(Hermes::Algebra::SimpleVector<double>* src, Hermes::Algebra::SimpleVector<double>* target, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  target->zero();
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    for(int i = 0; i < al_coarse.cnt; i++)
      target->set(al_fine.dof[i+1],  src->get(al_coarse.dof[i]));
  }
}

Hermes::Algebra::SimpleVector<double>* cut_off_ders(double* src_vector, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  SimpleVector<double>* vector = new SimpleVector<double>(space_coarse->get_num_dofs());
  vector->zero();
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    vector->set(al_coarse.dof[0],  src_vector[al_fine.dof[0]]);
  }
  return vector;
}

void add_ders(Hermes::Algebra::SimpleVector<double>* src, Hermes::Algebra::SimpleVector<double>* target, SpaceSharedPtr<double> space_coarse, SpaceSharedPtr<double> space_fine)
{
  target->zero();
  Element *e;
  for_all_active_elements(e, space_coarse->get_mesh())
  {
    AsmList<double> al_coarse, al_fine;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);

    target->set(al_fine.dof[0],  src->get(al_coarse.dof[0]));
  }
}
