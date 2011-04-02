#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Essential boundary conditions */

class EssentialBCNonConst : public EssentialBoundaryCondition {
public:
  EssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const 
         { return EssentialBoundaryCondition::BC_FUNCTION; }

  virtual scalar value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const {
    return -cos(M_PI*x);
  }
};

/* Weak forms */

class CustomWeakFormGeneral : public WeakForm
{
public:
  CustomWeakFormGeneral() : WeakForm(1)
  {
    add_matrix_form(new MatrixFormVolGeneral(0, 0));
    add_vector_form(new VectorFormVolGeneral(0));
    add_vector_form_surf(new VectorFormSurfGeneral(0, BDY_VERTICAL));
  }

  ~CustomWeakFormGeneral() {}

private:
  class MatrixFormVolGeneral : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolGeneral(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_SYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i=0; i < n; i++) {
        double x = e->x[i];
        double y = e->y[i];
        result += (a_11(x, y)*u->dx[i]*v->dx[i] +
                   a_12(x, y)*u->dy[i]*v->dx[i] +
                   a_21(x, y)*u->dx[i]*v->dy[i] +
                   a_22(x, y)*u->dy[i]*v->dy[i] +
                   a_1(x, y)*u->dx[i]*v->val[i] +
                   a_2(x, y)*u->dy[i]*v->val[i] +
                   a_0(x, y)*u->val[i]*v->val[i]) * wt[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
      // and test function plus two
    }
  private:
    // Problem parameters.
    double a_11(double x, double y) const { if (y > 0) return 1 + x*x + y*y; else return 1;}
    double a_22(double x, double y) const { if (y > 0) return 1; else return 1 + x*x + y*y;}
    double a_12(double x, double y) const { return 1; }
    double a_21(double x, double y) const { return 1;}
    double a_1(double x, double y) const { return 0.0;}
    double a_2(double x, double y) const { return 0.0;}
    double a_0(double x, double y) const { return 0.0;}
  };

  class VectorFormVolGeneral : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolGeneral(int i) : WeakForm::VectorFormVol(i) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }
  private:
    double rhs(double x, double y) const { return 1 + x*x + y*y;}
  };

  class VectorFormSurfGeneral : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfGeneral(int i, std::string area = HERMES_ANY) : WeakForm::VectorFormSurf(i, area) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (g_N(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }
  private:
    double g_N(double x, double y) const { return 0;}
  };
};
