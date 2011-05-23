/* Weak forms */

class WeakFormHeatTransferNewtonTimedep : public WeakForm<double>
{
public:
  WeakFormHeatTransferNewtonTimedep(double alpha, double tau, Solution<double>* sln_prev_time) 
          : WeakForm(1) 
  { 
    add_matrix_form(new MatrixFormVolHeatTransfer(0, 0, alpha, tau));
    add_vector_form(new VectorFormVolHeatTransfer(0, alpha, tau));
  };

private:
  class MatrixFormVolHeatTransfer : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, double alpha, double tau) 
          : MatrixFormVol<double>(i, j), alpha(alpha), tau(tau) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result -= wt[i] * (dlam_du<Real>(u_prev_newton->val[i]) * u->val[i] *
                           (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                           + lam<Real>(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext){
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext){
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const { 
      return 1 + pow(u, alpha); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) const { 
      return alpha*pow(u, alpha - 1); 
    }
    
    MatrixFormVol<double>* clone() {
      return new MatrixFormVolHeatTransfer(*this);
    }

    double alpha;
    double tau;
  };

  class VectorFormVolHeatTransfer : public VectorFormVol<double>
  {
  public:
    VectorFormVolHeatTransfer(int i, double alpha, double tau) 
          : VectorFormVol<double>(i), alpha(alpha), tau(tau) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result -= wt[i] * (lam<Real>(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] 
                           + u_prev_newton->dy[i] * v->dy[i])
		           - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext){
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext){
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) const {
      return 1.0;
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const { 
      return 1 + pow(u, alpha); 
    }

    VectorFormVol<double>* clone() {
      return new VectorFormVolHeatTransfer(*this);
    }

    double alpha;
    double tau;
  };
};

/* Essential boundary conditions */

class EssentialBCNonConst : public EssentialBoundaryCondition<double> {
public:
  EssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  //virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  virtual double value(double x, double y) const
  {
    return (x+10)*(y+10)/100.;
  }
};

/* Initial condition */
class InitialSolutionHeatTransfer : public ExactSolutionScalar<double>
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100.;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100.;
  }
};
