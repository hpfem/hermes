#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormHeatTransferPicard : public WeakForm
{
public:
  CustomWeakFormHeatTransferPicard(Solution* prev_iter_sln, double heat_src) : WeakForm(1) {
    MatrixFormVolHeatTransfer* matrix_form = new MatrixFormVolHeatTransfer(0, 0);
    matrix_form->ext.push_back(prev_iter_sln);
    add_matrix_form(matrix_form);

    add_vector_form(new DefaultVectorFormConst(0, heat_src));
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (lam<Real>(u_prev->val[i]) * (u->dx[i] * v->dx[i] 
                           + u->dy[i] * v->dy[i]));

      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // For any u, this function has to be positive.
    template<typename Real> 
    Real lam(Real u) const { 
      return 1 + pow(u, 4); 
    }
  };
};
