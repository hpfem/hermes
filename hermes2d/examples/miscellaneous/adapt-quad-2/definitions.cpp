class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(std::string material_1, double eps_1, std::string material_2, double eps_2, bool adapt_eval, int adapt_order_increase, double adapt_rel_error_tol) 
    : WeakForm(1) {
    MatrixFormLaplace* first_form = new MatrixFormLaplace(0, 0, material_1, eps_1);
    MatrixFormLaplace* second_form = new MatrixFormLaplace(0, 0, material_2, eps_2);
    if(adapt_eval) {
      first_form->adapt_order_increase = adapt_order_increase;
      second_form->adapt_order_increase = adapt_order_increase;
      first_form->adapt_rel_error_tol = adapt_rel_error_tol;
      second_form->adapt_rel_error_tol = adapt_rel_error_tol;
    }
    add_matrix_form(first_form);
    add_matrix_form(second_form);
  };

private:
  class MatrixFormLaplace : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormLaplace(int i, int j, std::string area, double a = 1) : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), a(a) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return a * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Constant multiplicating the Laplacian.
    const double a;
  };
};