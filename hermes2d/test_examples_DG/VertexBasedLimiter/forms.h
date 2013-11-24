#ifndef __VERTEX_BASED_FORMS
#define __VERTEX_BASED_FORMS

typedef double (*scalar_product_with_advection_direction)(double x, double y, double vx, double vy);
extern scalar_product_with_advection_direction advection_term;
extern bool only_x_der;

#pragma region Time derivative forms

class CustomMatrixFormVol : public MatrixFormVol<double>   
{
public:

  CustomMatrixFormVol(int i, int j, double factor) : MatrixFormVol<double>(i, j), factor(factor)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * v->val[i] * u->val[i];
    return result * this->factor;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
  {
    return v->val[0] * u->val[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVol(*this);
  }
  double factor;
};

class CustomVectorFormVol : public VectorFormVol<double>   
{
public:

  CustomVectorFormVol(int i, int prev, double factor) : VectorFormVol<double>(i), prev(prev), factor(factor)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * v->val[i] * ext[this->prev]->val[i];
    return result * this->factor;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
  {
    return v->val[0] * ext[this->prev]->val[0];
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVol(*this);
  }
  int prev;
  double factor;
};

#pragma endregion

#pragma region Error form

class ErrorFormSurf : public VectorFormSurf<double>
{
public:
  ErrorFormSurf(std::string area_to_set) : VectorFormSurf<double>(0)
  {
    this->set_area(area_to_set);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
      result += wt[i] * (ext[0]->val[i] - ext[1]->val[i]) * (ext[0]->val[i] - ext[1]->val[i]);

    return result;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return Ord(20);
  }

  VectorFormSurf<double>* clone() const
  {
    return new ErrorFormSurf(*this);
  }
};

class ErrorFormVol : public VectorFormVol<double>
{
public:
  ErrorFormVol() : VectorFormVol<double>(0)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
      result += wt[i] * (ext[0]->val[i] - ext[1]->val[i]) * (ext[0]->val[i] - ext[1]->val[i]);

    return result;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return Ord(20);
  }

  VectorFormVol<double>* clone() const
  {
    return new ErrorFormVol(*this);
  }
};

#pragma endregion

#pragma region Convection forms

class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  CustomMatrixFormVolConvection(int i, int j) : MatrixFormVol<double>(i, j) {}


  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * u->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
    return -result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return u->val[0] * v->dx[0] * e->nx[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVolConvection(*this);
  }
};

class CustomVectorFormVolConvection : public VectorFormVol<double>   
{
public:

  CustomVectorFormVolConvection(int i, int ext_i, double multiplier = 1.) : VectorFormVol<double>(i), ext_i(ext_i), multiplier(multiplier)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * ext[ext_i]->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
    return multiplier * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return ext[ext_i]->val[0] * v->dx[0] * e->nx[0];
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVolConvection(*this);
  }
  int ext_i;
  double multiplier;
};

class CustomMatrixFormInterfaceConvection : public MatrixFormDG<double>
{
public:
  CustomMatrixFormInterfaceConvection(int i, int j, bool local, bool inverse = false) : MatrixFormDG<double>(i, j), local(local), inverse(inverse)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    if(local)
    {
      if(inverse)
      {
        if(u->fn_central != NULL && v->fn_central != NULL)
          return 0.;
        if(u->fn_central == NULL && v->fn_central == NULL)
          return 0.;
      }
      else
      {
        if(u->fn_central == NULL && v->fn_central != NULL)
          return 0.;
        if(u->fn_central != NULL && v->fn_central == NULL)
          return 0.;
      }
    }
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

      double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);

      if(a_dot_n < 0)
        result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? u->val_neighbor[i] : 0.) * jump_v;
      if(a_dot_n >= 0)
        result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? 0. : u->val[i]) * jump_v;
    }

    return result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
  }

  MatrixFormDG<double>* clone() const
  {
    return new CustomMatrixFormInterfaceConvection(*this);
  }
  bool local, inverse;
};

class CustomVectorFormInterfaceConvection : public VectorFormDG<double>
{
public:
  CustomVectorFormInterfaceConvection(int i, int ext_i, bool on_K_in, bool on_K_out, double multiplier = -1.) : VectorFormDG<double>(i), ext_i(ext_i), on_K_in(on_K_in), on_K_out(on_K_out), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);
      if(a_dot_n > 0 && on_K_out)
        result += wt[i] * a_dot_n * ext[ext_i]->val[i] * v->val[i];
      if(a_dot_n < 0 && on_K_in)
        result += wt[i] * a_dot_n * ext[ext_i]->val_neighbor[i] * v->val[i];
    }
    return multiplier * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return Ord(3) * ext[ext_i]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterfaceConvection(*this);
  }

  int ext_i;
  bool on_K_in;
  bool on_K_out;
  double multiplier;
};

class CustomMatrixFormSurfConvection : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurfConvection(int i, int j)
    : MatrixFormSurf<double>(i, j)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double>* u, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);
      if(a_dot_n >= 0)
        result += wt[i] * u->val[i] * v->val[i] * a_dot_n;
    }
    return result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord>* u, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return u->val[0] * v->val[0] * e->x[0];
  }

  MatrixFormSurf<double>* clone() const
  {
    return new CustomMatrixFormSurfConvection(*this);
  }

};

class CustomVectorFormSurfConvection : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurfConvection(int i, int ext_bnd, bool on_K_in, bool on_K_out, double multiplier = -1.) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), on_K_in(on_K_in), on_K_out(on_K_out), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

      if(a_dot_n >= 0 && on_K_out)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
      if(a_dot_n < 0 && on_K_in)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
    }

    return multiplier * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return ext[this->ext_bnd]->val[0] * v->val[0] * e->nx[0];
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurfConvection(*this);
  }

  int ext_bnd;
  bool on_K_in;
  bool on_K_out;
  double multiplier;
};

#pragma endregion

#pragma region Diffusion forms

class CustomMatrixFormVolDiffusion : public MatrixFormVol<double>   
{
public:
  CustomMatrixFormVolDiffusion(int i, int j, double diffusivity) : MatrixFormVol<double>(i, j), diffusivity(diffusivity) {}

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->dx[i] * v->dx[i]);
    if(!only_x_der)
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dy[i] * v->dy[i]);

    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVolDiffusion(*this);
  }

  double diffusivity;
};

class CustomVectorFormVolDiffusion : public VectorFormVol<double>   
{
public:

  CustomVectorFormVolDiffusion(int i, int ext_i, double diffusivity, double multiplier = -1.) : VectorFormVol<double>(i), ext_i(ext_i), diffusivity(diffusivity), multiplier(multiplier)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext[ext_i]->dx[i] * v->dx[i]);
    if(!only_x_der)
      for (int i = 0; i < n; i++)
        result += wt[i] * (ext[ext_i]->dy[i] * v->dy[i]);
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return ext[ext_i]->val[0] * v->dx[0] * e->x[0] * Ord(3);
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVolDiffusion(*this);
  }
  int ext_i;
  double diffusivity;
  double multiplier;
};

class CustomMatrixFormInterfaceDiffusion : public MatrixFormDG<double>
{
public:
  CustomMatrixFormInterfaceDiffusion(int i, int j, bool local, double diffusivity, double s, double sigma, bool inverse = false) : MatrixFormDG<double>(i, j), local(local), inverse(inverse), diffusivity(diffusivity), s(s), sigma(sigma)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    if(local)
    {
      if(inverse)
      {
        if(u->fn_central != NULL && v->fn_central != NULL)
          return 0.;
        if(u->fn_central == NULL && v->fn_central == NULL)
          return 0.;
      }
      else
      {
        if(u->fn_central == NULL && v->fn_central != NULL)
          return 0.;
        if(u->fn_central != NULL && v->fn_central == NULL)
          return 0.;
      }
    }

    double result = 0.;

    switch(u->fn_central != NULL)
    {
    case true:
      switch(v->fn_central != NULL)
      {
      case true:
        // 0, 0
        for (int i = 0; i < n; i++) 
        {
          double jump_u = u->val[i];
          double jump_v = v->val[i];

          result -= 0.5 * wt[i] * (u->dx[i] * e->nx[i]) * jump_v;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (u->dy[i] * e->ny[i]) * jump_v;

          result -= 0.5 * wt[i] * (v->dx[i] * e->nx[i]) * jump_u * s;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (v->dy[i] * e->ny[i]) * jump_u * s;

          result += wt[i] * jump_v * jump_u * sigma;
        }
        break;
      case false:
        // 1, 0
        for (int i = 0; i < n; i++) 
        {
          double jump_u = u->val[i];
          double jump_v = -v->val_neighbor[i];

          result -= 0.5 * wt[i] * (u->dx[i] * e->nx[i]) * jump_v;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (u->dy[i] * e->ny[i]) * jump_v;

          result -= 0.5 * wt[i] * (v->dx_neighbor[i] * e->nx[i]) * jump_u * s;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (v->dy_neighbor[i] * e->ny[i]) * jump_u * s;

          result += wt[i] * jump_v * jump_u * sigma;
        }
      }
      break;
    case false:
      switch(v->fn_central != NULL)
      {
      case true:
        // 0, 1
        for (int i = 0; i < n; i++) 
        {
          double jump_u = -u->val_neighbor[i];
          double jump_v = v->val[i];

          result -= 0.5 * wt[i] * (u->dx_neighbor[i] * e->nx[i]) * jump_v;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (u->dy_neighbor[i] * e->ny[i]) * jump_v;

          result -= 0.5 * wt[i] * (v->dx[i] * e->nx[i]) * jump_u * s;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (v->dy[i] * e->ny[i]) * jump_u * s;

          result += wt[i] * jump_v * jump_u * sigma;
        }
        break;
      case false:
        // 1, 1
        for (int i = 0; i < n; i++) 
        {
          double jump_u = -u->val_neighbor[i];
          double jump_v = -v->val_neighbor[i];

          result -= 0.5 * wt[i] * (u->dx_neighbor[i] * e->nx[i]) * jump_v;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (u->dy_neighbor[i] * e->ny[i]) * jump_v;

          result -= 0.5 * wt[i] * (v->dx_neighbor[i] * e->nx[i]) * jump_u * s;
          if(!only_x_der)
            result -= 0.5 * wt[i] * (v->dy_neighbor[i] * e->ny[i]) * jump_u * s;

          result += wt[i] * jump_v * jump_u * sigma;
        }
      }
    }

    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
  }

  MatrixFormDG<double>* clone() const
  {
    return new CustomMatrixFormInterfaceDiffusion(*this);
  }
  bool local, inverse;
  double diffusivity;
  double s;
  double sigma;
};

class CustomVectorFormInterfaceDiffusion : public VectorFormDG<double>
{
public:
  CustomVectorFormInterfaceDiffusion(int i, int ext_i, double diffusivity, double s, double sigma, double multiplier = -1.) : VectorFormDG<double>(i), ext_i(ext_i), diffusivity(diffusivity), s(s), sigma(sigma), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = .5 * (ext[this->ext_i]->dx[i] + ext[this->ext_i]->dx_neighbor[i]);
      double dy = only_x_der ? 0. : .5 * (ext[this->ext_i]->dy[i] + ext[this->ext_i]->dy_neighbor[i]);
      result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      double dx_test = .5 * v->dx[i];
      double dy_test = only_x_der ? 0. : .5 * v->dy[i];
      double jump = ext[this->ext_i]->val[i] - ext[this->ext_i]->val_neighbor[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * jump * s;
      result += wt[i] * sigma * jump * v->val[i];
    }
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return ext[this->ext_i]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterfaceDiffusion(*this);
  }

  int ext_i;
  double s, sigma;
  double diffusivity;
  double multiplier;
};

class CustomVectorFormInterfaceDiffusionOffDiag : public VectorFormDG<double>
{
public:
  CustomVectorFormInterfaceDiffusionOffDiag(int i, int ext_i, double diffusivity, double s, double sigma, double multiplier = -1.) : VectorFormDG<double>(i), ext_i(ext_i), diffusivity(diffusivity), s(s), sigma(sigma), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = .5 * ext[this->ext_i]->dx_neighbor[i];
      double dy = only_x_der ? 0. : .5 * ext[this->ext_i]->dy_neighbor[i];
      result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      double dx_test = .5 * v->dx[i];
      double dy_test = only_x_der ? 0. : .5 * v->dy[i];
      double jump = - ext[this->ext_i]->val_neighbor[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * jump * s;
      result += wt[i] * sigma * jump * v->val[i];
    }
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return ext[this->ext_i]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterfaceDiffusionOffDiag(*this);
  }

  int ext_i;
  double s, sigma;
  double diffusivity;
  double multiplier;
};

class CustomMatrixFormSurfDiffusion : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurfDiffusion(int i, int j, double diffusivity, double s, double sigma, std::string inlet)
    : MatrixFormSurf<double>(i, j), diffusivity(diffusivity), s(s), sigma(sigma)
  {
    this->set_area(inlet);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double>* u, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = u->dx[i];
      double dy = only_x_der ? 0. : u->dy[i];
      result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      double dx_test = v->dx[i];
      double dy_test = only_x_der ? 0. : v->dy[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * u->val[i] * s;
      result += wt[i] * sigma * v->val[i] * u->val[i];
    }
    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord>* u, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return u->val[0] * v->val[0] * e->x[0] * e->nx[0] * Ord(3);
  }

  MatrixFormSurf<double>* clone() const
  {
    return new CustomMatrixFormSurfDiffusion(*this);
  }
  double diffusivity;
  double s, sigma;
};

class CustomVectorFormSurfDiffusion : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurfDiffusion(int i, int ext_bnd, double diffusivity, double s, double sigma, std::string inlet, bool add_grad_u = true, double multiplier = -1.) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), diffusivity(diffusivity), s(s), sigma(sigma), add_grad_u(add_grad_u), multiplier(multiplier)
  {
    this->set_area(inlet);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = ext[this->ext_bnd]->dx[i];
      double dy = only_x_der ? 0. : ext[this->ext_bnd]->dy[i];
      double dx_test = v->dx[i];
      double dy_test = only_x_der ? 0. : v->dy[i];
      if(add_grad_u)
        result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * ext[this->ext_bnd]->val[i] * s;
      result += wt[i] * sigma * v->val[i] * ext[this->ext_bnd]->val[i];
    }
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return ext[this->ext_bnd]->val[0] * v->val[0] * e->nx[0] * Ord(3);
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurfDiffusion(*this);
  }

  int ext_bnd;
  double diffusivity;
  double s, sigma, multiplier;
  bool add_grad_u;
};

class CustomVectorFormSurfDiffusionNeumann : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurfDiffusionNeumann(int i, int ext_bnd, std::string outlet, double diffusivity, MeshFunctionSharedPtr<double> exact_solution) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), diffusivity(diffusivity)
  {
    this->set_area(outlet);
    this->set_ext(exact_solution);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
      result += wt[i] * v->val[i] * (e->nx[i] * ext[this->ext_bnd]->dx[i] + e->ny[i] * ext[this->ext_bnd]->dy[i]);
    return diffusivity * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return Ord(20);
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurfDiffusionNeumann(*this);
  }

  int ext_bnd;
  double diffusivity;
};
#pragma endregion
#endif
