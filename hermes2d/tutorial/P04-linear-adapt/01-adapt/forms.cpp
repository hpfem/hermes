#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Basic electrostatic class - need improvement
// TODO: cannot handle BCs
class WeakFormElectrostatic : public WeakForm
{
public:
  // Problem parameters.
  double EPS0;

  // Relative electric permittivity
  std::map<int, double> eps_r;
  // Charge density
  std::map<int, double> rho;

  WeakFormElectrostatic() : WeakForm(1)
  {
    EPS0 = 8.863e-12;

    add_matrix_form(new MatrixFormVolElectrostatic(0, 0));
    add_vector_form(new VectorFormVolElectrostatic(0));
  }

  class MatrixFormVolElectrostatic : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolElectrostatic(int i, int j) : WeakForm::MatrixFormVol(i, j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return static_cast<WeakFormElectrostatic *>(wf)->eps_r[e->elem_marker]
          * static_cast<WeakFormElectrostatic *>(wf)->EPS0
          * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolElectrostatic : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolElectrostatic(int i) : WeakForm::VectorFormVol(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return static_cast<WeakFormElectrostatic *>(wf)->rho[e->elem_marker] * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };
};


// An example of ose of Electrostatic class
class WeakFormTutorial : public WeakFormElectrostatic
{
public:
  WeakFormTutorial() : WeakFormElectrostatic()
  {
    // Problem parameters.
    const int MATERIAL_1 = 1;
    const int MATERIAL_2 = 2;

    eps_r[MATERIAL_1] = 1.0;
    eps_r[MATERIAL_2] = 10.0;

    rho[MATERIAL_1] = 0.0;
    rho[MATERIAL_2] = 0.0;
  }
};
