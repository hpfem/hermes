class CustomExactFunction
{
public:
  CustomExactFunction(double t_init, double final_time) : t_init(t_init), final_time(final_time) {}
  
  double temp_ext(double t) {
    return t_init + 10. * sin(2 * M_PI * t / final_time);
  }

  // Members.
  double t_init;
  double final_time;
};

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double heatcap, double rho, double lambda, double tau, double alpha, Solution* prev_time_sln, CustomExactFunction* exact_temp_ext) : WeakForm(1) {
    MatrixFormVol* matrix_form_vol = new MatrixFormVol(0, 0, heatcap, rho, lambda, tau, alpha);
    VectorFormVol* vector_form_vol = new VectorFormVol(0, heatcap, rho, lambda, tau, alpha);
    vector_form_vol->ext.push_back(prev_time_sln);
    MatrixFormSurf* matrix_form_surf = new MatrixFormSurf(0, 0, heatcap, rho, lambda, tau, alpha);
    VectorFormSurf* vector_form_surf = new VectorFormSurf(0, heatcap, rho, lambda, tau, alpha, exact_temp_ext);

    add_matrix_form(matrix_form_vol);
    add_vector_form(vector_form_vol);
    add_matrix_form_surf(matrix_form_surf);
    add_vector_form_surf(vector_form_surf);
  };

private:
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j, double heatcap, double rho, double lambda, double tau, double alpha) : WeakForm::MatrixFormVol(i, j), heatcap(heatcap),
                  rho(rho), lambda(lambda), tau(tau), alpha(alpha) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return heatcap * rho * int_u_v<Real, Scalar>(n, wt, u, v) / tau +
         lambda * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double heatcap, rho, lambda, tau, alpha;
  };

  class VectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVol(int i, double heatcap, double rho, double lambda, double tau, double alpha) : WeakForm::VectorFormVol(i), heatcap(heatcap),
                  rho(rho), lambda(lambda), tau(tau), alpha(alpha) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      return heatcap * rho * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / tau;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    double heatcap, rho, lambda, tau, alpha;
  };

  class MatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormSurf(int i, int j, double heatcap, double rho, double lambda, double tau, double alpha) : WeakForm::MatrixFormSurf(i, j), heatcap(heatcap),
                  rho(rho), lambda(lambda), tau(tau), alpha(alpha) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return lambda * alpha * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
    
    // Members.
    double heatcap, rho, lambda, tau, alpha;
  };

  class VectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurf(int i, double heatcap, double rho, double lambda, double tau, double alpha, CustomExactFunction* exact_temp_ext) : WeakForm::VectorFormSurf(i), heatcap(heatcap),
                  rho(rho), lambda(lambda), tau(tau), alpha(alpha), exact_temp_ext(exact_temp_ext) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      return lambda * alpha * exact_temp_ext->temp_ext(wf->get_current_time()) * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Time-dependent exterior temperature.
    CustomExactFunction* exact_temp_ext;

    // Members.
    double heatcap, rho, lambda, tau, alpha;
  };

};