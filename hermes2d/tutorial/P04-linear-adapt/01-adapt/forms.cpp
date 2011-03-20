#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Basic electrostatic class - need improvement
class WeakFormElectrostatic : public WeakForm
{
public:
  // Problem parameters.
  double EPS0;

  // Relative electric permittivity
  std::map<std::string, double> eps_r;
  // Charge density
  std::map<std::string, double> rho;

  WeakFormElectrostatic() : WeakForm(1)
  {
    EPS0 = 8.863e-12;

    add_matrix_form(new MatrixFormVolElectrostatic(0, 0));
    add_vector_form(new VectorFormVolElectrostatic(0));
  }

  class MatrixFormVolElectrostatic : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolElectrostatic(int i, int j) : WeakForm::MatrixFormVol(i, j) {
    adapt_eval = false;
    sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return static_cast<WeakFormElectrostatic *>(wf)->eps_r[wf->get_element_markers_conversion()->get_user_marker(e->elem_marker)]
          * static_cast<WeakFormElectrostatic *>(wf)->EPS0
          * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }
  };

  class VectorFormVolElectrostatic : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolElectrostatic(int i) : WeakForm::VectorFormVol(i) {
    adapt_eval = false;
    }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return static_cast<WeakFormElectrostatic *>(wf)->rho[wf->get_element_markers_conversion()->get_user_marker(e->elem_marker)] * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return int_v<Ord, Ord>(n, wt, v);
    }
  };
};


// An example of use of Electrostatic class
class WeakFormElectrostaticTutorial : public WeakFormElectrostatic
{
public:
  WeakFormElectrostaticTutorial() : WeakFormElectrostatic()
  {
    // Set material properties    
    // Element markers.
    eps_r["1"] = 1.0;
    eps_r["2"] = 10.0;

    rho["1"] = 0.0;
    rho["2"] = 0.0;
  }
};
