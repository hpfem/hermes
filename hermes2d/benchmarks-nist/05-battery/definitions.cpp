#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string omega_1, std::string omega_2, 
                        std::string omega_3, std::string omega_4, 
                        std::string omega_5, std::string bdy_left, 
                        std::string bdy_top, std::string bdy_right, 
                        std::string bdy_bottom) : WeakForm(1),
  
  omega_1(omega_1), omega_2(omega_2), omega_3(omega_3), 
  omega_4(omega_4), omega_5(omega_5), 

  p_1(25.0),
  p_2(7.0),
  p_3(5.0),
  p_4(0.2),
  p_5(0.05),

  q_1(25.0),
  q_2(0.8),
  q_3(0.0001),
  q_4(0.2),
  q_5(0.05),

  f_1(0.0),
  f_2(1.0),
  f_3(1.0),
  f_4(0.0),
  f_5(0.0),

  bdy_left(bdy_left), 
  bdy_top(bdy_top), 
  bdy_right(bdy_right), 
  bdy_bottom(bdy_bottom),

  c_left(0.0),
  c_top(1.0),
  c_right(2.0),
  c_bottom(3.0),

  g_n_left(0.0),
  g_n_top(3.0),
  g_n_right(2.0),
  g_n_bottom(1.0)

  {
    add_matrix_form(new CustomMatrixFormVol(0, 0));
    add_vector_form(new CustomVectorFormVol(0));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_bottom));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_right));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_top));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_bottom));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_top));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_left));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_right));
  };

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_SYM) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      double p = 0, q = 0;
      // integration order calculation.
      if(e->elem_marker == -9999)
        p = q = 1;
      else {
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_1) {
            p = static_cast<CustomWeakFormPoisson*>(wf)->p_1;
            q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
        }
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_2) {
            p = static_cast<CustomWeakFormPoisson*>(wf)->p_2;
            q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
        }
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_3) {
            p = static_cast<CustomWeakFormPoisson*>(wf)->p_3;
            q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
        }
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_4) {
            p = static_cast<CustomWeakFormPoisson*>(wf)->p_4;
            q = static_cast<CustomWeakFormPoisson*>(wf)->q_4;
        }
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_5) {
            p = static_cast<CustomWeakFormPoisson*>(wf)->p_5;
            q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
        }
      }
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (p * u->dx[i] * v->dx[i] + q * u->dy[i] * v->dy[i]);
      return result;
    };

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      double f=0;
      if(e->elem_marker == -9999)
        f = 1;
      else {
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_1)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_1;
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_2)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_2;
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_3)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_3;
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_4)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_4;
        if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) 
           == static_cast<CustomWeakFormPoisson*>(wf)->omega_5)
            f = static_cast<CustomWeakFormPoisson*>(wf)->f_5;
      }

      return f * int_v<Real, Scalar>(n, wt, v);
    };

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class CustomMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    CustomMatrixFormSurf(int i, int j, std::string marker) 
          : WeakForm::MatrixFormSurf(i, j, marker) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        Real x = e->x[i];
        Real y = e->y[i];
        Scalar p = 0.0;
        Scalar q = 0.0;
        Scalar c = 1.0;
        if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_left) {
          if (x == 0.0) {
            if ((y >= 0.0 && y <= 0.8)||(y >= 23.2 && y <= 24.0)) {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
            }
            if ((y >= 1.6 && y <= 3.6)||(y >= 18.8 && y <= 21.2)) {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_2; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_2;
            }
            if (y >= 3.6 && y <= 18.8) {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_3; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_3;
            }
            if ((y >= 0.8 && y <= 1.6)||(y >= 21.2 && y <= 23.2)) {
              p = static_cast<CustomWeakFormPoisson*>(wf)->p_5; 
              q = static_cast<CustomWeakFormPoisson*>(wf)->q_5;
            }
          }
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_left;
        }
        if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_right) {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_right;
        }

        if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_bottom) {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_bottom;
        }

        if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_top) {
          p = static_cast<CustomWeakFormPoisson*>(wf)->p_1; 
          q = static_cast<CustomWeakFormPoisson*>(wf)->q_1;
          c = static_cast<CustomWeakFormPoisson*>(wf)->c_top;
        }
        result += wt[i] * (p * u->dx[i] * v->val[i] - q * u->dy[i] * v->val[i] 
                  + c * u->val[i] * v->val[i]);
      }
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return wt[0] * (u->dx[0] * v->val[0] - u->dy[0] * v->val[0] +  u->val[0] * v->val[0]);
    }
  };

  class CustomVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurf(int i, std::string marker) 
          : WeakForm::VectorFormSurf(i, marker) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Scalar g = 1.0;
      if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_left) {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_left; 
      }
      if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_right) {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_right;
      }

      if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_bottom) {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_bottom;
      }

      if(this->area == static_cast<CustomWeakFormPoisson*>(wf)->bdy_top) {
        g = static_cast<CustomWeakFormPoisson*>(wf)->g_n_top;
      }
      return g * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  // Members.
  const std::string omega_1;
  const std::string omega_2;
  const std::string omega_3;
  const std::string omega_4;
  const std::string omega_5;

  const double p_1;
  const double p_2;
  const double p_3;
  const double p_4;
  const double p_5;

  const double q_1;
  const double q_2;
  const double q_3;
  const double q_4;
  const double q_5;

  const double f_1;
  const double f_2;
  const double f_3;
  const double f_4;
  const double f_5;

  // Boundary markers.
  const std::string bdy_left;
  const std::string bdy_top;
  const std::string bdy_right;
  const std::string bdy_bottom;

  // Boundary condition coefficients for the four sides.
  const double c_left;
  const double c_top;
  const double c_right;
  const double c_bottom;

  const double g_n_left;
  const double g_n_top;
  const double g_n_right;
  const double g_n_bottom;
};

