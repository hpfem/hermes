// The first part of the file dontains forms for the Newton's
// method. Forms for Picard are in the second part.

/*** INITIAL CONDITION ***/

class InitialSolutionRichards : public ExactSolutionScalar
{
public:
  InitialSolutionRichards(Mesh* mesh, double constant) 
         : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0.0;
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }

  // Value.
  double constant;
};

class ExactSolutionPoisson : public ExactSolutionScalar
{
public:
  ExactSolutionPoisson(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return x*x +y*y;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 2*x;
    dy = 2*y;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return x*x +y*y;
  }
};

/*** NEWTON ***/

#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormRichardsNewtonEuler : public WeakForm
{
public:
  WeakFormRichardsNewtonEuler(ConstitutiveRelations* relations, double tau, Solution* prev_time_sln) 
               : WeakForm(1) {
    JacobianFormNewtonEuler* jac_form = new JacobianFormNewtonEuler(0, 0, relations, tau);
    jac_form->ext.push_back(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonEuler* res_form = new ResidualFormNewtonEuler(0, relations, tau);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonEuler : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormNewtonEuler(int i, int j, ConstitutiveRelations* relations, double tau) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];

      for (int i = 0; i < n; i++)
        result += wt[i] * (
        relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                          );
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormNewtonEuler : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormNewtonEuler(int i, ConstitutiveRelations* relations, double tau) 
               : WeakForm::VectorFormVol(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           -relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          );
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };
};


class WeakFormRichardsNewtonCrankNicolson : public WeakForm
{
public:
  WeakFormRichardsNewtonCrankNicolson(ConstitutiveRelations* relations, double tau, Solution* prev_time_sln) : WeakForm(1) {
    JacobianFormNewtonCrankNicolson* jac_form = new JacobianFormNewtonCrankNicolson(0, 0, relations, tau);
    jac_form->ext.push_back(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonCrankNicolson* res_form = new ResidualFormNewtonCrankNicolson(0, relations, tau);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonCrankNicolson : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormNewtonCrankNicolson(int i, int j, ConstitutiveRelations* relations, double tau) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * 0.5 * ( // implicit Euler part:
		             relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                           )
                + wt[i] * 0.5 * ( // explicit Euler part, 
		             relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                           );
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormNewtonCrankNicolson : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormNewtonCrankNicolson(int i, ConstitutiveRelations* relations, double tau) : WeakForm::VectorFormVol(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      double result = 0;
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * ( // implicit Euler part
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          )
                + wt[i] * 0.5 * ( // explicit Euler part
		           relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_time->val[i], atoi(elem_marker.c_str())) * h_prev_time->dy[i] * v->val[i]
		           );
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };
};


class WeakFormRichardsPicardEuler : public WeakForm
{
public:
  WeakFormRichardsPicardEuler(ConstitutiveRelations* relations, double tau, Solution* prev_picard_sln, Solution* prev_time_sln) : WeakForm(1) {
    JacobianFormPicardEuler* jac_form = new JacobianFormPicardEuler(0, 0, relations, tau);
    jac_form->ext.push_back(prev_picard_sln);
    add_matrix_form(jac_form);

    ResidualFormPicardEuler* res_form = new ResidualFormPicardEuler(0, relations, tau);
    res_form->ext.push_back(prev_picard_sln);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormPicardEuler : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormPicardEuler(int i, int j, ConstitutiveRelations* relations, double tau) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      double result = 0;
      Func<double>* h_prev_picard = ext->fn[0];

      for (int i = 0; i < n; i++) {
        result += wt[i] * (  relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                             + relations->K(h_prev_picard->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             - relations->dKdh(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]);
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormPicardEuler : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormPicardEuler(int i, ConstitutiveRelations* relations, double tau) : WeakForm::VectorFormVol(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = wf->get_element_markers_conversion()->get_user_marker(e->elem_marker);
      double result = 0;
      Func<double>* h_prev_picard = ext->fn[0];
      Func<double>* h_prev_time = ext->fn[1];
      for (int i = 0; i < n; i++) 
        result += wt[i] * relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * h_prev_time->val[i] * v->val[i] / tau;
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };
};

class RichardsEssentialBC : public EssentialBoundaryCondition {
public:

  RichardsEssentialBC(std::string marker, double h_elevation, double pulse_end_time, double h_init, double startup_time) :
  EssentialBoundaryCondition(Hermes::vector<std::string>()), h_elevation(h_elevation), pulse_end_time(pulse_end_time), h_init(h_init), startup_time(startup_time)
  {
    markers.push_back(marker);
  }

  ~RichardsEssentialBC() {}

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    if (current_time < startup_time)
      return h_init + current_time/startup_time*(h_elevation-h_init);
    else if (current_time > pulse_end_time)
      return h_init;
    else
      return h_elevation;
  }

  // Member.
  double h_elevation;
  double pulse_end_time;
  double h_init;
  double startup_time;
};

/*** EXTRAS ***/

//Debugging matrix printer.
bool printmatrix(double** A, int n, int m){
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++){
      printf(" %lf ", A[i][j]) ;
    }
    printf(" \n");
  }
  printf("----------------------------------\n");
  return true;
}


//Debugging vector printer.
bool printvector(double* vect, int n){
  for (int i=0; i<n; i++){
    printf(" %lf ", vect[i]);
  }
  printf("\n");
  printf("----------------------------------\n");
  return true;
}
