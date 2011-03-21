#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakForm : public WeakForm
{
public:

  CustomWeakForm(std::string left_bottom_bnd_part) : WeakForm(1) {
    add_matrix_form(new MatrixFormVol(0, 0));
    add_vector_form(new VectorFormVol(0));
    add_matrix_form_surf(new MatrixFormSurface(0, 0));
    add_matrix_form_surf(new MatrixFormInterface(0, 0));
    add_vector_form_surf(new VectorFormSurface(0, left_bottom_bnd_part));


  };

private:
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += -wt[i] * u->val[i] * static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(e->x[i], e->y[i], v->dx[i], v->dy[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVol(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * F(e->x[i], e->y[i]) * v->val[i];
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Real>
    Real F(Real x, Real y) {
      return 0;
    }
  };

  class MatrixFormSurface : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormSurface(int i, int j) : WeakForm::MatrixFormSurf(i, j, H2D_DG_BOUNDARY_EDGE) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
  
      for (int i = 0; i < n; i++) {
        Real x = e->x[i], y = e->y[i];
        Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(x, y, e->nx[i], e->ny[i]);
        result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], 0, a_dot_n) * v->val[i];
      }
  
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class MatrixFormInterface : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormInterface(int i, int j) : WeakForm::MatrixFormSurf(i, j, H2D_DG_INNER_EDGE) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
  
      for (int i = 0; i < n; i++) {
        Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(e->x[i], e->y[i], e->nx[i], e->ny[i]);
        Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
        result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
      }
  
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormSurface : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurface(int i, std::string left_bottom_bnd_part) : WeakForm::VectorFormSurf(i, H2D_DG_BOUNDARY_EDGE), left_bottom_bnd_part(left_bottom_bnd_part) { }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      double result = 0;
      for (int i = 0; i < n; i++) {
        double x = e->x[i], y = e->y[i];
        double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(x, y, e->nx[i], e->ny[i]);
        // Function values for Dirichlet boundary conditions.
        result += -wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0, g<double,double>(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker), x, y), a_dot_n) * v->val[i];
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      Ord result = 0;
      for (int i = 0; i < n; i++)
        result += -wt[i] * v->val[i];
      return result;
    }

    template<typename Real>
    Real F(Real x, Real y) {
      return 0;
    }

    template<typename Real, typename Scalar>
    Scalar g(std::string ess_bdy_marker, Real x, Real y) {
      if (ess_bdy_marker == left_bottom_bnd_part) return 1; else return 0;
    }
    
    // Member.
    std::string left_bottom_bnd_part;
  };
  
  double calculate_a_dot_v(double x, double y, double vx, double vy) {
    double norm = std::max<double>(1e-12, sqrt(sqr(x) + sqr(y)));
    return -y/norm*vx + x/norm*vy;
  }

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) {
    return Ord(10);
  }

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) {
    return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib); 
  }

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) {
    return a_dot_n * (u_cent + u_neib); 
  }
};