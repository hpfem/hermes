#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

class WeakFormPoissonNIST05 : public WeakForm
{
public:
  WeakFormPoissonNIST05(double parameter) : WeakForm(1)
  {
    add_matrix_form(new MatrixFormVolPoisson(0, 0));
    
    VectorFormVolPoisson* wfp= new VectorFormVolPoisson(0);
    wfp->parameter = parameter;
    add_vector_form(wfp);
  };

private:
  class MatrixFormVolPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPoisson(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      double p, q;
      // integration order calculation.
      if(e->elem_marker = -8888)
        p = q = 1;
      else
        switch(wf->element_markers_conversion->get_user_marker(e->elem_marker)) {
          case OMEGA_1:
            p = P_1;
            q = Q_1;
            break;
          case OMEGA_2:
            p = P_2;
            q = Q_2;
            break;
          case OMEGA_3:
            p = P_3;
            q = Q_3;
            break;
          case OMEGA_4:
            p = P_4;
            q = Q_4;
            break;
          case OMEGA_5:
            p = P_5;
            q = Q_5;
            break;
        }
    };

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    const std::string OMEGA_1;
    const std::string OMEGA_2;
    const std::string OMEGA_3;
    const std::string OMEGA_4;
    const std::string OMEGA_5;

    const double P_1;
    const double P_2;
    const double P_3;
    const double P_4;
    const double P_5;

    const double Q_1;
    const double Q_2;
    const double Q_3;
    const double Q_4;
    const double Q_5;
  };

  class VectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPoisson(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      double p;
      if(e->elem_marker = -8888)
        p = 1;
      else
        switch(wf->element_markers_conversion->get_user_marker(e->elem_marker)) {
          case OMEGA_1:
            p = F_1;
            break;
          case OMEGA_2:
            p = F_2;
            break;
          case OMEGA_3:
            p = F_3;
            break;
          case OMEGA_4:
            p = F_4;
            break;
          case OMEGA_5:
            p = F_5;
            break;
        }
      return p * int_v<Real, Scalar>(n, wt, v);
    };

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
     
    // Members.
    const std::string OMEGA_1;
    const std::string OMEGA_2;
    const std::string OMEGA_3;
    const std::string OMEGA_4;
    const std::string OMEGA_5;

    const double F_1;
    const double F_2;
    const double F_3;
    const double F_4;
    const double F_5;
  };

  // Members.
  const std::string OMEGA_1 = "1";
  const std::string OMEGA_2 = "2";
  const std::string OMEGA_3 = "3";
  const std::string OMEGA_4 = "4";
  const std::string OMEGA_5 = "5";

  const double P_1 = 25.0;
  const double P_2 = 7.0;
  const double P_3 = 5.0;
  const double P_4 = 0.2;
  const double P_5 = 0.05;

  const double Q_1 = 25.0;
  const double Q_2 = 0.8;
  const double Q_3 = 0.0001;
  const double Q_4 = 0.2;
  const double Q_5 = 0.05;

  const double F_1 = 0.0;
  const double F_2 = 1.0;
  const double F_3 = 1.0;
  const double F_4 = 0.0;
  const double F_5 = 0.0;

  // Boundary markers.
  const std::string BDY_LEFT = "1";
  const std::string BDY_TOP = "2";
  const std::string BDY_RIGHT = "3";
  const std::string BDY_BOTTOM = "4";

  // Boundary condition coefficients for the four sides.
  const double C_LEFT = 0.0;
  const double C_TOP = 1.0;
  const double C_RIGHT = 2.0;
  const double C_BOTTOM = 3.0;

  const double G_N_LEFT = 0.0;
  const double G_N_TOP = 3.0;
  const double G_N_RIGHT = 2.0;
  const double G_N_BOTTOM = 1.0;
};

class DirichletFunctionBoundaryConditionExact : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryConditionExact(std::string marker) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()) 
  {
    markers.push_back(marker);
  };
  
  ~DirichletFunctionBoundaryConditionExact() {};

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return 0; //exact_solution->fn(x, y);
  };
};

































// Weak forms
template<typename Real, typename Scalar>
Scalar bilinear_form_surf_left(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                               Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double x = e->x[i];
    double y = e->y[i];
    double P = 0.0;
    double Q = 0.0;

    if (x == 0.0)
    {
      if ((y >= 0.0 && y <= 0.8)||(y >= 23.2 && y <= 24.0))
      {
        P = P_1; 
        Q = Q_1;
      }
      if ((y >= 1.6 && y <= 3.6)||(y >= 18.8 && y <= 21.2))
      {
        P = P_2; 
        Q = Q_2;
      }
      if (y >= 3.6 && y <= 18.8)
      {
        P = P_3; 
        Q = Q_3;
      }
      if ((y >= 0.8 && y <= 1.6)||(y >= 21.2 && y <= 23.2))
      {
        P = P_5; 
        Q = Q_5;
      }
    }
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_LEFT * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_right(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                                Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_RIGHT * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                              Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_TOP * u->val[i] * v->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v,
                                 Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    double P = 25.0;
    double Q = 25.0;
    result += wt[i] * (P * u->dx[i] * v->val[i] - Q * u->dy[i] * v->val[i] + C_BOTTOM * u->val[i] * v->val[i]);
  }
  return result;
}

Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
//return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_left(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_LEFT * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_right(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_RIGHT * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_TOP * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return G_N_BOTTOM * int_v<Real, Scalar>(n, wt, v);
}
