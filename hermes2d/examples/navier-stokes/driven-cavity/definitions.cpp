#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::SurfaceMatrixForms;
using namespace WeakFormsH1::SurfaceVectorForms;

/* Weak forms */

class WeakFormDrivenCavity : public WeakForm
{
public:
  WeakFormDrivenCavity(double Re, std::string bdy_top, double time_step, Solution* x_vel_previous_time, 
                         Solution* y_vel_previous_time) 
    : WeakForm(3), Re(Re), time_step(time_step), x_vel_previous_time(x_vel_previous_time), 
                y_vel_previous_time(y_vel_previous_time) {
    /* Jacobian terms - first velocity equation */
    // Time derivative in the first velocity equation.
    add_matrix_form(new DefaultLinearMass(0, 0, 1./time_step)); 
    // Laplacian divided by Re in the first velocity equation.
    add_matrix_form(new DefaultLinearDiffusion(0, 0, 1./Re));
    // First part of the convective term in the first velocity equation.
    BilinearFormNonsymVel_0_0* nonsym_vel_form_0_0 = new BilinearFormNonsymVel_0_0(0, 0);
    add_matrix_form(nonsym_vel_form_0_0);
    // Second part of the convective term in the first velocity equation.
    BilinearFormNonsymVel_0_1* nonsym_vel_form_0_1 = new BilinearFormNonsymVel_0_1(0, 1);
    add_matrix_form(nonsym_vel_form_0_1);
    // Pressure term in the first velocity equation.
    BilinearFormNonsymXVelPressure* nonsym_velx_pressure_form = new BilinearFormNonsymXVelPressure(0, 2);
    add_matrix_form(nonsym_velx_pressure_form);

    /* Jacobian terms - second velocity equation, continuity equation */
    // Time derivative in the second velocity equation.
    add_matrix_form(new DefaultLinearMass(1, 1, 1./time_step));
    // Laplacian divided by Re in the second velocity equation.
    add_matrix_form(new DefaultLinearDiffusion(1, 1, 1./Re));
    // First part of the convective term in the second velocity equation.
    BilinearFormNonsymVel_1_0* nonsym_vel_form_1_0 = new BilinearFormNonsymVel_1_0(1, 0);
    add_matrix_form(nonsym_vel_form_1_0);
    // Second part of the convective term in the second velocity equation.
    BilinearFormNonsymVel_1_1* nonsym_vel_form_1_1 = new BilinearFormNonsymVel_1_1(1, 1);
    add_matrix_form(nonsym_vel_form_1_1);
    // Pressure term in the second velocity equation.
    BilinearFormNonsymYVelPressure* nonsym_vely_pressure_form = new BilinearFormNonsymYVelPressure(1, 2);
    add_matrix_form(nonsym_vely_pressure_form);

    /* Residual - volumetric */
    // First velocity equation.
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Re, time_step);
    F_0->ext = Hermes::vector<MeshFunction*>(x_vel_previous_time);
    add_vector_form(F_0);
    // Second velocity equation.
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Re, time_step);
    F_1->ext = Hermes::vector<MeshFunction*>(y_vel_previous_time);
    add_vector_form(F_1);
    // Continuity equation. 
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);
  };

  class BilinearFormNonsymVel_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      Func<scalar>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] 
			   + u->val[i] * xvel_prev_newton->dx[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
    Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] 
			   + u->val[i] * xvel_prev_newton->dx[i] * v->val[i]);
      return result;
    }
  };

  class BilinearFormNonsymVel_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * xvel_prev_newton->dy[i] * v->val[i] ;

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;

      Func<Ord>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
 
      return result;
    }
  };


  class BilinearFormNonsymVel_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* yvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* yvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      return result;
    }
  };


  class BilinearFormNonsymVel_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      Func<scalar>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] 
                           + u->val[i] * yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] 
                           + u->val[i] * yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }
  };

  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    // The antisym flag is used here to generate a term in the continuity equation.
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANTISYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return - int_u_dvdx<double, scalar>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
  };


  class BilinearFormNonsymYVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    // The antisym flag is used here to generate a term in the continuity equation.
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANTISYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return - int_u_dvdy<double, scalar>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }
  };

  class VectorFormNS_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_0(int i, double Re, double time_step) : WeakForm::VectorFormVol(i), Re(Re), time_step(time_step) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_time = ext->fn[0];  
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  
      Func<scalar>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Re 
                          - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step 
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_time = ext->fn[0];  
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Re 
                          - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step 
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    double Re, time_step;
  };

  class VectorFormNS_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_1(int i, double Re, double time_step) 
      : WeakForm::VectorFormVol(i), Re(Re), time_step(time_step) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* yvel_prev_time = ext->fn[0];  
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  
      Func<scalar>* p_prev_newton = u_ext[2];
      Func<scalar>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Re 
                          - (p_prev_newton->val[i] * v->dy[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] 
			      + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* yvel_prev_time = ext->fn[0];  
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Re 
                  - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    double Re, time_step;
  };

  class VectorFormNS_2 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_2(int i) : WeakForm::VectorFormVol(i) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] + yvel_prev_newton->dy[i]) * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] + yvel_prev_newton->dy[i]) * v->val[i];
      return result;
    }
  };

  class VectorFormNS_3 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_3(int i, double time_step) : WeakForm::VectorFormVol(i), time_step(time_step) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* temp_prev_time = ext->fn[0];  
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  
      Func<scalar>* temp_prev_newton = u_ext[3];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (((temp_prev_newton->val[i] - temp_prev_time->val[i]) / time_step
                           + xvel_prev_newton->val[i] * temp_prev_newton->dx[i] 
			   + yvel_prev_newton->val[i] * temp_prev_newton->dy[i]) * v->val[i]
	                   + temp_prev_newton->dx[i] * v->dx[i] + temp_prev_newton->dy[i] * v->dy[i]
			   );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* temp_prev_time = ext->fn[0];  
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* temp_prev_newton = u_ext[3];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (((temp_prev_newton->val[i] - temp_prev_time->val[i]) / time_step
                           + xvel_prev_newton->val[i] * temp_prev_newton->dx[i] 
			   + yvel_prev_newton->val[i] * temp_prev_newton->dy[i]) * v->val[i]
	                   + temp_prev_newton->dx[i] * v->dx[i] + temp_prev_newton->dy[i] * v->dy[i]
			   );
      return result;
    }

    private:
      double time_step;
  };

  class BilinearFormNonsymTemp_3_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymTemp_3_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dx[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dx[i] * v->val[i];
      return result;
    }
  };

  class BilinearFormNonsymTemp_3_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymTemp_3_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dy[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dy[i] * v->val[i];
      return result;
    }
  };

  class BilinearFormNonsymTemp_3_3 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymTemp_3_3(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      Func<scalar>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->val[i] * u->dx[i] 
                           + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->val[i] * u->dx[i] 
                           + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i];
      return result;
    }
  };

  class CustomResidualSurfConst : public WeakForm::VectorFormSurf
  {
  public:
    CustomResidualSurfConst(int i, scalar coeff = 1.0,
                             GeomType gt = HERMES_PLANAR)
           : WeakForm::VectorFormSurf(i), coeff(coeff), gt(gt) { }
    CustomResidualSurfConst(int i, std::string area, scalar coeff = 1.0,
                             GeomType gt = HERMES_PLANAR)
           : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * u_ext[3]->val[i] * v->val[i];
      }
      return coeff * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::VectorFormSurf* clone() {
      return new CustomResidualSurfConst(*this);
    }

    private:
      scalar coeff;
      GeomType gt;
  };


protected:
  double Re, time_step;
  Solution* x_vel_previous_time;
  Solution* y_vel_previous_time;
};

