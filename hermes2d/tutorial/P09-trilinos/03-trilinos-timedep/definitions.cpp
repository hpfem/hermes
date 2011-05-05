class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(Hermes::vector<std::string> newton_boundaries, double heatcap, 
    double rho, double tau, double lambda, double alpha, double temp_ext, Solution* sln_prev_time, bool JFNK = false) : WeakForm(1, JFNK)
  {
    // Jacobian forms - volumetric.
    add_matrix_form(new JacobianFormVol(0, 0, heatcap, rho, lambda, tau));

    // Jacobian forms - surface.
    add_matrix_form_surf(new JacobianFormSurf(0, 0, newton_boundaries, alpha, lambda));


    // Residual forms - volumetric.
    ResidualFormVol* res_form = new ResidualFormVol(0, heatcap, rho, lambda, tau);
    res_form->ext.push_back(sln_prev_time);
    add_vector_form(res_form);

    // Residual forms - surface.
    add_vector_form_surf(new ResidualFormSurf(0, newton_boundaries, alpha, lambda, temp_ext));
  }

  ~CustomWeakForm() {}

private:
  class JacobianFormVol : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormVol(int i, int j, double heatcap, double rho, double lambda, double tau) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM),
    heatcap(heatcap), rho(rho), lambda(lambda), tau(tau) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (heatcap * rho * u->val[i] * v->val[i] / tau
                         + lambda * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double heatcap, rho, lambda, tau;
  };

  class JacobianFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    JacobianFormSurf(int i, int j, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda) : WeakForm::MatrixFormSurf(i, j, newton_boundaries),
      alpha(alpha), lambda(lambda) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * alpha * lambda * u->val[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double alpha, lambda;
  };

  class ResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormVol(int i, double heatcap, double rho, double lambda, double tau) : WeakForm::VectorFormVol(i),
    heatcap(heatcap), rho(rho), lambda(lambda), tau(tau)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (heatcap * rho * (u_ext[0]->val[i] - ext->fn[0]->val[i]) * v->val[i] / tau
		           + lambda * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double heatcap, rho, lambda, tau;
  };
  
  class ResidualFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    ResidualFormSurf(int i, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda, double temp_ext) : WeakForm::VectorFormSurf(i, newton_boundaries),
    alpha(alpha), lambda(lambda), temp_ext(temp_ext)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * alpha * lambda * (u_ext[0]->val[i] - temp_ext) * v->val[i];
      return result;  
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double alpha, lambda, temp_ext;
  };
};