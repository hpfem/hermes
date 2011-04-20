#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "../euler_util.h"

class EulerEquationsWeakFormExplicit : public WeakForm
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicit(NumericalFlux* num_flux, double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
  std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
  Solution* prev_density, Solution* prev_density_vel_x, Solution* prev_density_vel_y, Solution* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
  WeakForm(num_of_equations), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
  energy_ext(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)), euler_fluxes(new EulerFluxes(kappa)) {
    add_matrix_form(new EulerEquationsBilinearFormTime(0));
    add_matrix_form(new EulerEquationsBilinearFormTime(1));
    add_matrix_form(new EulerEquationsBilinearFormTime(2));
    add_matrix_form(new EulerEquationsBilinearFormTime(3));
    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    if(!fvm_only) {
      add_vector_form(new EulerEquationsLinearFormDensity());
      add_vector_form(new EulerEquationsLinearFormDensityVelX(kappa));
      add_vector_form(new EulerEquationsLinearFormDensityVelY(kappa));
      add_vector_form(new EulerEquationsLinearFormEnergy(kappa));
    }

    add_vector_form_surf(new EulerEquationsLinearFormInterface(0, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(1, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(2, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(3, num_flux));

    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_bottom_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_bottom_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_bottom_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_bottom_marker, num_flux));

    if(solid_wall_bottom_marker != solid_wall_top_marker) {
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_top_marker, num_flux));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_top_marker, num_flux));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_top_marker, num_flux));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_top_marker, num_flux));
    }
    else
      warning("Are you sure that solid wall top and bottom markers should coincide?");

    add_vector_form_surf(new EulerEquationsLinearFormInlet(0, inlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(1, inlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(2, inlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(3, inlet_marker, num_flux));

    add_vector_form_surf(new EulerEquationsLinearFormOutlet(0, outlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(1, outlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(2, outlet_marker, num_flux));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(3, outlet_marker, num_flux));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol.size();vector_form_i++) {
      vfvol.at(vector_form_i)->ext.push_back(prev_density);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfvol.at(vector_form_i)->ext.push_back(prev_energy);
    }

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf.size();vector_form_i++) {
      vfsurf.at(vector_form_i)->ext.push_back(prev_density);
      vfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfsurf.at(vector_form_i)->ext.push_back(prev_energy);
    }
  };

  void set_time_step(double tau) {
    this->tau = tau;
  }

  double get_tau() const {
    return tau;
  }

  // Destructor.
  ~EulerEquationsWeakFormExplicit() {};
protected:
  class EulerEquationsBilinearFormTime : public WeakForm::MatrixFormVol
  {
  public:
    EulerEquationsBilinearFormTime(int i) : WeakForm::MatrixFormVol(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class EulerEquationsLinearFormDensity : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensity() : WeakForm::VectorFormVol(0) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, 
                       ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
         ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
      ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class EulerEquationsLinearFormDensityVelX : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensityVelX(double kappa)
     : WeakForm::VectorFormVol(1), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    double kappa;
  };

  class EulerEquationsLinearFormDensityVelY : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensityVelY(double kappa) 
         : WeakForm::VectorFormVol(2), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  
    double kappa;
  };

  class EulerEquationsLinearFormEnergy : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormEnergy(double kappa) 
         : WeakForm::VectorFormVol(3), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                        * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                        * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                        * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                        * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
                        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                        * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    double kappa;
  };

  class EulerEquationsLinearFormTime : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormTime(int i) 
         : WeakForm::VectorFormVol(i), component_i(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      return int_u_v<Real, Scalar>(n, wt, ext->fn[component_i], v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    int component_i;
  };

  class EulerEquationsLinearFormInterface : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormInterface(int i, NumericalFlux* num_flux) 
         : WeakForm::VectorFormSurf(i, H2D_DG_INNER_EDGE), element(i), num_flux(num_flux) {}

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      double result = 0;
      double w_L[4], w_R[4];
      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->get_val_central(i);
        w_R[0] = ext->fn[0]->get_val_neighbor(i);
    
        w_L[1] = ext->fn[1]->get_val_central(i);
        w_R[1] = ext->fn[1]->get_val_neighbor(i);

        w_L[2] = ext->fn[2]->get_val_central(i);
        w_R[2] = ext->fn[2]->get_val_neighbor(i);

        w_L[3] = ext->fn[3]->get_val_central(i);
        w_R[3] = ext->fn[3]->get_val_neighbor(i);

        result -= wt[i] * v->val[i] 
                        * num_flux->numerical_flux_i(element, w_L, w_R, e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormSolidWall : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormSolidWall(int i, std::string marker, NumericalFlux* num_flux) 
         : WeakForm::VectorFormSurf(i, marker), element(i), num_flux(num_flux) {}

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      double result = 0;
      for (int i = 0;i < n;i++) {
        double w_L[4];
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        result -= wt[i] * v->val[i] * num_flux->numerical_flux_solid_wall_i(element, 
                  w_L, e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormInlet : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormInlet(int i, std::string marker, NumericalFlux* num_flux) 
         : WeakForm::VectorFormSurf(i, marker), element(i), num_flux(num_flux) {}

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                 ExtData<scalar> *ext) const {
      double result = 0;
      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) {
        // Left (inner) state from the previous time level solution.
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];
    
        w_B[0] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext;
        w_B[1] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v1_ext;
        w_B[2] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v2_ext;
        w_B[3] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->energy_ext;
        
        result -= wt[i] * v->val[i] * num_flux->numerical_flux_inlet_i(element, 
                  w_L, w_B, e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormOutlet : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormOutlet(int i, std::string marker, NumericalFlux* num_flux) : 
    WeakForm::VectorFormSurf(i, marker), element(i), num_flux(num_flux) {}

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      double result = 0;
      double w_L[4];
      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];
        result -= wt[i] * v->val[i] 
                        * num_flux->numerical_flux_outlet_i(element, w_L, 
                          static_cast<EulerEquationsWeakFormExplicit*>(wf)->pressure_ext, 
                          e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };
  // Members.
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double pressure_ext;
  double energy_ext;
  double tau;

  EulerFluxes* euler_fluxes;

  friend class EulerEquationsWeakFormExplicitCoupled;
  friend class EulerEquationsWeakFormExplicitMultiComponent;
  friend class EulerEquationsWeakFormSemiImplicitMultiComponent;
  friend class EulerEquationsWeakFormSemiImplicitCoupled;
};

class EulerEquationsWeakFormExplicitMultiComponent : public WeakForm
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicitMultiComponent(NumericalFlux* num_flux, double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
        std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
        Solution* prev_density, Solution* prev_density_vel_x, Solution* prev_density_vel_y, Solution* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
        WeakForm(num_of_equations), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
        energy_ext(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)), euler_fluxes(new EulerFluxes(kappa))
  {
    Hermes::vector<std::pair<unsigned int, unsigned int> > matrix_coordinates;
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(0, 0));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(1, 1));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(2, 2));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(3, 3));

    Hermes::vector<unsigned int> vector_coordinates;
    vector_coordinates.push_back(0);
    vector_coordinates.push_back(1);
    vector_coordinates.push_back(2);
    vector_coordinates.push_back(3);

    add_multicomponent_matrix_form(new EulerEquationsBilinearFormTime(matrix_coordinates));

    if(!fvm_only)
      add_multicomponent_vector_form(new EulerEquationsLinearForm(vector_coordinates, kappa));

    add_multicomponent_vector_form(new EulerEquationsLinearFormTime(vector_coordinates));
    
    add_multicomponent_vector_form_surf(new EulerEquationsLinearFormInterface(vector_coordinates, num_flux));

    add_multicomponent_vector_form_surf(new EulerEquationsLinearFormSolidWall(vector_coordinates, 
                                        solid_wall_bottom_marker, num_flux));
    
    if(solid_wall_bottom_marker != solid_wall_top_marker)
      add_multicomponent_vector_form_surf(new EulerEquationsLinearFormSolidWall(vector_coordinates, 
                                          solid_wall_top_marker, num_flux));

    add_multicomponent_vector_form_surf(new EulerEquationsLinearFormInlet(vector_coordinates, 
                                        inlet_marker, num_flux));

    add_multicomponent_vector_form_surf(new EulerEquationsLinearFormOutlet(vector_coordinates, 
                                        outlet_marker, num_flux));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol_mc.size();vector_form_i++) {
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_energy);
    }

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf_mc.size();vector_form_i++) {
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_energy);
    }
  };

  void set_time_step(double tau) {
    this->tau = tau;
  }

  double get_tau() const {
    return tau;
  }

  // Destructor.
  ~EulerEquationsWeakFormExplicitMultiComponent() {};
protected:
  class EulerEquationsBilinearFormTime : public WeakForm::MultiComponentMatrixFormVol
  {
  public:
    EulerEquationsBilinearFormTime(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates) : WeakForm::MultiComponentMatrixFormVol(coordinates) {}

  void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_n = int_u_v<double, double>(n, wt, u, v);
      
      result.push_back(result_n);
      result.push_back(result_n);
      result.push_back(result_n);
      result.push_back(result_n);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return int_u_v<Ord, Ord>(n, wt, u, v);
    }
  };

  class EulerEquationsLinearForm : public WeakForm::MultiComponentVectorFormVol
  {
  public:
    EulerEquationsLinearForm(Hermes::vector<unsigned int> coordinates, double kappa) 
                             : WeakForm::MultiComponentVectorFormVol(coordinates), kappa(kappa) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;
      for (int i = 0;i < n;i++) {
        result_0 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        
        result_1 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        
        result_2 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
     
        result_3 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dx[i];
        result_3 += wt[i] * ext->fn[0]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dy[i];
        result_3 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dx[i];
        result_3 += wt[i] * ext->fn[1]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_3 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_3 += wt[i] * ext->fn[2]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dy[i];
        result_3 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_3 += wt[i] * ext->fn[3]->val[i] 
                          * (static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
      }
      result.push_back(result_0 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_1 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_2 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_3 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0] * v->val[0] * v->val[0] * v->val[0];
    }

    double kappa;
  };
  
  class EulerEquationsLinearFormTime : public WeakForm::MultiComponentVectorFormVol
  {
  public:
    EulerEquationsLinearFormTime(Hermes::vector<unsigned int> coordinates) 
         : WeakForm::MultiComponentVectorFormVol(coordinates) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[0], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[1], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[2], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[3], v));
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return int_u_v<Ord, Ord>(n, wt, v, v);
    }
  };

  class EulerEquationsLinearFormInterface : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    EulerEquationsLinearFormInterface(Hermes::vector<unsigned int> coordinates, 
                                      NumericalFlux* num_flux) 
         : WeakForm::MultiComponentVectorFormSurf(coordinates, H2D_DG_INNER_EDGE), num_flux(num_flux) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;
      double w_L[4], w_R[4];
      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->get_val_central(i);
        w_R[0] = ext->fn[0]->get_val_neighbor(i);
    
        w_L[1] = ext->fn[1]->get_val_central(i);
        w_R[1] = ext->fn[1]->get_val_neighbor(i);

        w_L[2] = ext->fn[2]->get_val_central(i);
        w_R[2] = ext->fn[2]->get_val_neighbor(i);

        w_L[3] = ext->fn[3]->get_val_central(i);
        w_R[3] = ext->fn[3]->get_val_neighbor(i);

        double flux[4];
        num_flux->numerical_flux(flux, w_L, w_R, e->nx[i], e->ny[i]);

#ifdef H2D_EULER_NUM_FLUX_TESTING
        double flux_testing_num_flux[4];
        double flux_testing_num_flux_conservativity_1[4];
        double flux_testing_num_flux_conservativity_2[4];
        double flux_testing_flux_1[4];
        double flux_testing_flux_2[4];
        double flux_testing_flux[4];

        num_flux->numerical_flux(flux_testing_num_flux, w_L, w_L, e->nx[i], e->ny[i]);

        num_flux->numerical_flux(flux_testing_num_flux_conservativity_1, w_L, w_R, e->nx[i], e->ny[i]);
        num_flux->numerical_flux(flux_testing_num_flux_conservativity_2, w_R, w_L, -e->nx[i], -e->ny[i]);

       for(unsigned int flux_i = 0;flux_i < 4;flux_i++)
          if(std::abs(flux_testing_num_flux_conservativity_1[flux_i] + flux_testing_num_flux_conservativity_2[flux_i]) > 1E-4)
            if(std::abs((flux_testing_num_flux_conservativity_1[flux_i] + flux_testing_num_flux_conservativity_2[flux_i]) / flux_testing_num_flux_conservativity_1[flux_i]) > 1E-6)
              info("Flux is not conservative.");

        //num_flux->Q(w_L, w_L, e->nx[i], e->ny[i]);
        
        EulerEquationsLinearForm form(Hermes::vector<unsigned int>(), num_flux->kappa);
        flux_testing_flux_1[0] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_1[1] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_1[2] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_1[3] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_2[0] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_2[1] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_2[2] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];

        flux_testing_flux_2[3] = form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_0<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[0]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_1<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[1]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_2<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[2]
          + form.(static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_3<double>(w_L[0], w_L[1], w_L[2], w_L[3]) * w_L[3];


          flux_testing_flux[0] = flux_testing_flux_1[0] * e->nx[i] + flux_testing_flux_2[0] * e->ny[i];
          flux_testing_flux[1] = flux_testing_flux_1[1] * e->nx[i] + flux_testing_flux_2[1] * e->ny[i];
          flux_testing_flux[2] = flux_testing_flux_1[2] * e->nx[i] + flux_testing_flux_2[2] * e->ny[i];
          flux_testing_flux[3] = flux_testing_flux_1[3] * e->nx[i] + flux_testing_flux_2[3] * e->ny[i];

        //num_flux->Q_inv(flux_testing_flux_1, flux_testing_flux_1, e->nx[i], e->ny[i]);

        for(unsigned int flux_i = 0;flux_i < 4;flux_i++)
          if(std::abs(flux_testing_num_flux[flux_i] - flux_testing_flux[flux_i]) > 1E-8)
            if(std::abs((flux_testing_num_flux[flux_i] - flux_testing_flux[flux_i]) / flux_testing_num_flux[flux_i]) > 1E-6)
              info("Flux is not consistent.");
#endif

        result_0 -= wt[i] * v->val[i] * flux[0];
        result_1 -= wt[i] * v->val[i] * flux[1];
        result_2 -= wt[i] * v->val[i] * flux[2];
        result_3 -= wt[i] * v->val[i] * flux[3];
      }
      result.push_back(result_0 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_1 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_2 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_3 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormSolidWall : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    EulerEquationsLinearFormSolidWall(Hermes::vector<unsigned int> coordinates, 
                                      std::string marker, NumericalFlux* num_flux) 
         : WeakForm::MultiComponentVectorFormSurf(coordinates, marker), num_flux(num_flux) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;
      for (int i = 0;i < n;i++) {
        double w_L[4];
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        double flux[4];
        num_flux->numerical_flux_solid_wall(flux, w_L, e->nx[i], e->ny[i]);

        result_0 -= wt[i] * v->val[i] * flux[0];
        result_1 -= wt[i] * v->val[i] * flux[1];
        result_2 -= wt[i] * v->val[i] * flux[2];
        result_3 -= wt[i] * v->val[i] * flux[3];
      }
      result.push_back(result_0 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_1 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_2 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_3 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormInlet : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    EulerEquationsLinearFormInlet(Hermes::vector<unsigned int> coordinates, 
                                  std::string marker, NumericalFlux* num_flux) : 
        WeakForm::MultiComponentVectorFormSurf(coordinates, marker), num_flux(num_flux) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;
      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) {
        // Left (inner) state from the previous time level solution.
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];
    
        w_B[0] = static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->rho_ext;
        w_B[1] = static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->v1_ext;
        w_B[2] = static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->v2_ext;
        w_B[3] = static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->energy_ext;
    
        double flux[4];
        num_flux->numerical_flux(flux, w_L, w_B, e->nx[i], e->ny[i]);

        result_0 -= wt[i] * v->val[i] * flux[0];
        result_1 -= wt[i] * v->val[i] * flux[1];
        result_2 -= wt[i] * v->val[i] * flux[2];
        result_3 -= wt[i] * v->val[i] * flux[3];
      }
      result.push_back(result_0 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_1 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_2 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_3 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormOutlet : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    EulerEquationsLinearFormOutlet(Hermes::vector<unsigned int> coordinates, 
                                   std::string marker, NumericalFlux* num_flux) : 
    WeakForm::MultiComponentVectorFormSurf(coordinates, marker), num_flux(num_flux) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;
      double w_L[4];
      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        double flux[4];
        num_flux->numerical_flux_outlet(flux, w_L, 
          static_cast<EulerEquationsWeakFormExplicitMultiComponent*>(wf)->pressure_ext, e->nx[i], e->ny[i]);
        result_0 -= wt[i] * v->val[i] * flux[0];
        result_1 -= wt[i] * v->val[i] * flux[1];
        result_2 -= wt[i] * v->val[i] * flux[2];
        result_3 -= wt[i] * v->val[i] * flux[3];
      }

      result.push_back(result_0 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_1 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_2 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
      result.push_back(result_3 * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0];
    }

    // Members.
    NumericalFlux* num_flux;
  };

  // Members.
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double pressure_ext;
  double energy_ext;
  double tau;
  EulerFluxes* euler_fluxes;
};

class EulerEquationsWeakFormSemiImplicitMultiComponent : public WeakForm
{
public:
  // Constructor.
  EulerEquationsWeakFormSemiImplicitMultiComponent(NumericalFlux* num_flux, double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
        std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
        Solution* prev_density, Solution* prev_density_vel_x, Solution* prev_density_vel_y, Solution* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
        WeakForm(num_of_equations), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
        energy_ext(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)), euler_fluxes(new EulerFluxes(kappa))
  {
    Hermes::vector<std::pair<unsigned int, unsigned int> > matrix_coordinates;
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(0, 0));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(1, 1));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(2, 2));
    matrix_coordinates.push_back(std::pair<unsigned int, unsigned int>(3, 3));

    Hermes::vector<std::pair<unsigned int, unsigned int> > matrix_coordinates_full;
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(0, 0));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(0, 1));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(0, 2));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(0, 3));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(1, 0));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(1, 1));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(1, 2));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(1, 3));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(2, 0));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(2, 1));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(2, 2));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(2, 3));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(3, 0));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(3, 1));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(3, 2));
    matrix_coordinates_full.push_back(std::pair<unsigned int, unsigned int>(3, 3));

    Hermes::vector<unsigned int> vector_coordinates;
    vector_coordinates.push_back(0);
    vector_coordinates.push_back(1);
    vector_coordinates.push_back(2);
    vector_coordinates.push_back(3);

    add_multicomponent_matrix_form(new EulerEquationsBilinearFormTime(matrix_coordinates));

    if(!fvm_only)
      add_multicomponent_matrix_form(new EulerEquationsBilinearForm(matrix_coordinates_full, kappa));
    
    add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSurfSemiImplicit(matrix_coordinates_full, kappa));
    add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSurfSemiImplicit2(matrix_coordinates_full, kappa));

    add_multicomponent_vector_form(new EulerEquationsLinearFormTime(vector_coordinates));
    
    add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(matrix_coordinates_full, 
                                        solid_wall_bottom_marker, kappa));
    
    if(solid_wall_bottom_marker != solid_wall_top_marker)
      add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(matrix_coordinates_full, 
                                          solid_wall_top_marker, kappa));

    add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInlet(matrix_coordinates_full, 
                                        inlet_marker, kappa));
    add_multicomponent_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitOutlet(matrix_coordinates_full, 
                                        outlet_marker, kappa));

    add_multicomponent_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInlet(vector_coordinates, 
                                        inlet_marker, kappa));
    //add_multicomponent_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(vector_coordinates, 
    //                                    outlet_marker, kappa));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol_mc.size();vector_form_i++) {
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfvol_mc.at(vector_form_i)->ext.push_back(prev_energy);
    }

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf_mc.size();vector_form_i++) {
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfsurf_mc.at(vector_form_i)->ext.push_back(prev_energy);
    }
    
    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfvol_mc.size();matrix_form_i++) {
      mfvol_mc.at(matrix_form_i)->ext.push_back(prev_density);
      mfvol_mc.at(matrix_form_i)->ext.push_back(prev_density_vel_x);
      mfvol_mc.at(matrix_form_i)->ext.push_back(prev_density_vel_y);
      mfvol_mc.at(matrix_form_i)->ext.push_back(prev_energy);
    }

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfsurf_mc.size();matrix_form_i++) {
      mfsurf_mc.at(matrix_form_i)->ext.push_back(prev_density);
      mfsurf_mc.at(matrix_form_i)->ext.push_back(prev_density_vel_x);
      mfsurf_mc.at(matrix_form_i)->ext.push_back(prev_density_vel_y);
      mfsurf_mc.at(matrix_form_i)->ext.push_back(prev_energy);
    }
  };

  void set_time_step(double tau) {
    this->tau = tau;
  }

  double get_tau() const {
    return tau;
  }

  // Destructor.
  ~EulerEquationsWeakFormSemiImplicitMultiComponent() {};
protected:
  class EulerEquationsBilinearFormTime : public WeakForm::MultiComponentMatrixFormVol
  {
  public:
    EulerEquationsBilinearFormTime(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates) : WeakForm::MultiComponentMatrixFormVol(coordinates) {}

  void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_n = int_u_v<double, double>(n, wt, u, v);
      
      result.push_back(result_n);
      result.push_back(result_n);
      result.push_back(result_n);
      result.push_back(result_n);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return int_u_v<Ord, Ord>(n, wt, u, v);
    }
  };

  class EulerEquationsBilinearForm : public WeakForm::MultiComponentMatrixFormVol
  {
  public:
    EulerEquationsBilinearForm(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double kappa) 
                             : WeakForm::MultiComponentMatrixFormVol(coordinates), kappa(kappa) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      for (int i = 0;i < n;i++) {
        result_0_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_0_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_0_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_0_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_0_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        
        result_1_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_1_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_1_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_1_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_1_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        
        result_2_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_2_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_2_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_2_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_2_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
     
        result_3_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dx[i];
        result_3_0 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_0<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dy[i];
        result_3_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dx[i];
        result_3_1 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_1<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
        result_3_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_3_2 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_2<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
                          * v->dy[i];
        result_3_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_1_3_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dx[i];
        result_3_3 += wt[i] * u->val[i]
                          * (static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf))->euler_fluxes->A_2_3_3<double>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) 
                          * v->dy[i];
      }
      result.push_back(result_0_0 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_0_1 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_0_2 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_0_3 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());

      result.push_back(result_1_0 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_1_1 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_1_2 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_1_3 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());

      result.push_back(result_2_0 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_2_1 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_2_2 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_2_3 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());

      result.push_back(result_3_0 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_3_1 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_3_2 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(result_3_3 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return u->val[0] * u->val[0] + v->val[0] * v->val[0];
    }

    double kappa;
  };

  class EulerEquationsMatrixFormSurfSemiImplicit : public WeakForm::MultiComponentMatrixFormSurf
  {
  public:
    EulerEquationsMatrixFormSurfSemiImplicit(Hermes::vector<std::pair<unsigned int, 
                                                unsigned int> >coordinates, double kappa) 
    : WeakForm::MultiComponentMatrixFormSurf(coordinates, H2D_DG_INNER_EDGE), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext, 
                 Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      double result_L[4];
      double result_R[4];

      double w_L[4], w_R[4];

      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->get_val_central(i);
        w_R[0] = ext->fn[0]->get_val_neighbor(i);
    
        w_L[1] = ext->fn[1]->get_val_central(i);
        w_R[1] = ext->fn[1]->get_val_neighbor(i);

        w_L[2] = ext->fn[2]->get_val_central(i);
        w_R[2] = ext->fn[2]->get_val_neighbor(i);

        w_L[3] = ext->fn[3]->get_val_central(i);
        w_R[3] = ext->fn[3]->get_val_neighbor(i);

        double e_1_1_1[4] = {1, 0, 0, 0};
        double e_2_2_1[4] = {0, 1, 0, 0};
        double e_3_3_1[4] = {0, 0, 1, 0};
        double e_4_4_1[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];

        num_flux->P_plus(P_plus_1, w_L, e_1_1_1, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_2, w_L, e_2_2_1, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_3, w_L, e_3_3_1, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_4, w_L, e_4_4_1, e->nx[i], e->ny[i]);

        result_0_0 += wt[i] * (P_plus_1[0] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_1 += wt[i] * (P_plus_2[0] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_2 += wt[i] * (P_plus_3[0] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_3 += wt[i] * (P_plus_4[0] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_1_0 += wt[i] * (P_plus_1[1] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_1 += wt[i] * (P_plus_2[1] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_2 += wt[i] * (P_plus_3[1] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_3 += wt[i] * (P_plus_4[1] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_2_0 += wt[i] * (P_plus_1[2] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_1 += wt[i] * (P_plus_2[2] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_2 += wt[i] * (P_plus_3[2] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_3 += wt[i] * (P_plus_4[2] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_3_0 += wt[i] * (P_plus_1[3] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_1 += wt[i] * (P_plus_2[3] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_2 += wt[i] * (P_plus_3[3] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_3 += wt[i] * (P_plus_4[3] * u->get_val_central(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
      }

      result.push_back(result_0_0);
      result.push_back(result_0_1);
      result.push_back(result_0_2);
      result.push_back(result_0_3);

      result.push_back(result_1_0);
      result.push_back(result_1_1);
      result.push_back(result_1_2);
      result.push_back(result_1_3);

      result.push_back(result_2_0);
      result.push_back(result_2_1);
      result.push_back(result_2_2);
      result.push_back(result_2_3);

      result.push_back(result_3_0);
      result.push_back(result_3_1);
      result.push_back(result_3_2);
      result.push_back(result_3_3);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * u->get_val_central(i) * v->get_val_central(i);
        result += wt[i] * u->get_val_neighbor(i) * v->get_val_neighbor(i);
      }
      return result;
    }

    StegerWarmingNumericalFlux* num_flux;
  };
  
  class EulerEquationsMatrixFormSurfSemiImplicit2 : public WeakForm::MultiComponentMatrixFormSurf
  {
  public:
    EulerEquationsMatrixFormSurfSemiImplicit2(Hermes::vector<std::pair<unsigned int, 
                                                unsigned int> >coordinates, double kappa) 
    : WeakForm::MultiComponentMatrixFormSurf(coordinates, H2D_DG_INNER_EDGE), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext, 
                 Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      double result_L[4];
      double result_R[4];

      double w_L[4], w_R[4];

      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->get_val_central(i);
        w_R[0] = ext->fn[0]->get_val_neighbor(i);
    
        w_L[1] = ext->fn[1]->get_val_central(i);
        w_R[1] = ext->fn[1]->get_val_neighbor(i);

        w_L[2] = ext->fn[2]->get_val_central(i);
        w_R[2] = ext->fn[2]->get_val_neighbor(i);

        w_L[3] = ext->fn[3]->get_val_central(i);
        w_R[3] = ext->fn[3]->get_val_neighbor(i);
        
        double e_1_1_2[4] = {1, 0, 0, 0};
        double e_2_2_2[4] = {0, 1, 0, 0};
        double e_3_3_2[4] = {0, 0, 1, 0};
        double e_4_4_2[4] = {0, 0, 0, 1};
              
        double P_minus_1[4];
        double P_minus_2[4];
        double P_minus_3[4];
        double P_minus_4[4];
      
        num_flux->P_minus(P_minus_1, w_R, e_1_1_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_2, w_R, e_2_2_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_3, w_R, e_3_3_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_4, w_R, e_4_4_2, e->nx[i], e->ny[i]);

        result_0_0 += wt[i] * (P_minus_1[0] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_1 += wt[i] * (P_minus_2[0] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_2 += wt[i] * (P_minus_3[0] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_3 += wt[i] * (P_minus_4[0] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_1_0 += wt[i] * (P_minus_1[1] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_1 += wt[i] * (P_minus_2[1] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_2 += wt[i] * (P_minus_3[1] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_3 += wt[i] * (P_minus_4[1] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_2_0 += wt[i] * (P_minus_1[2] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_1 += wt[i] * (P_minus_2[2] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_2 += wt[i] * (P_minus_3[2] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_3 += wt[i] * (P_minus_4[2] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_3_0 += wt[i] * (P_minus_1[3] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_1 += wt[i] * (P_minus_2[3] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_2 += wt[i] * (P_minus_3[3] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_3 += wt[i] * (P_minus_4[3] * u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

      }
      
      result.push_back(result_0_0);
      result.push_back(result_0_1);
      result.push_back(result_0_2);
      result.push_back(result_0_3);

      result.push_back(result_1_0);
      result.push_back(result_1_1);
      result.push_back(result_1_2);
      result.push_back(result_1_3);

      result.push_back(result_2_0);
      result.push_back(result_2_1);
      result.push_back(result_2_2);
      result.push_back(result_2_3);

      result.push_back(result_3_0);
      result.push_back(result_3_1);
      result.push_back(result_3_2);
      result.push_back(result_3_3);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      for (int i = 0;i < n;i++) {
        result += wt[i] * u->get_val_central(i) * v->get_val_central(i);
        result += wt[i] * u->get_val_neighbor(i) * v->get_val_neighbor(i);
      }
      return result;
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsMatrixFormSemiImplicitInlet : public WeakForm::MultiComponentMatrixFormSurf
  {
  public:
    EulerEquationsMatrixFormSemiImplicitInlet(Hermes::vector<std::pair<unsigned int, 
                                                unsigned int> >coordinates, 
                                  std::string marker, double kappa) 
    : WeakForm::MultiComponentMatrixFormSurf(coordinates, marker), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext, 
               Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        double e_1_1[4] = {1, 0, 0, 0};
        double e_2_2[4] = {0, 1, 0, 0};
        double e_3_3[4] = {0, 0, 1, 0};
        double e_4_4[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];
      
        num_flux->P_plus(P_plus_1, w_L, e_1_1, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_2, w_L, e_2_2, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_3, w_L, e_3_3, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_4, w_L, e_4_4, e->nx[i], e->ny[i]);

        result_0_0 += wt[i] * P_plus_1[0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_1 += wt[i] * P_plus_2[0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_2 += wt[i] * P_plus_3[0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_3 += wt[i] * P_plus_4[0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_1_0 += wt[i] * P_plus_1[1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_1 += wt[i] * P_plus_2[1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_2 += wt[i] * P_plus_3[1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_3 += wt[i] * P_plus_4[1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_2_0 += wt[i] * P_plus_1[2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_1 += wt[i] * P_plus_2[2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_2 += wt[i] * P_plus_3[2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_3 += wt[i] * P_plus_4[2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_3_0 += wt[i] * P_plus_1[3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_1 += wt[i] * P_plus_2[3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_2 += wt[i] * P_plus_3[3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_3 += wt[i] * P_plus_4[3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
      }

      result.push_back(result_0_0);
      result.push_back(result_0_1);
      result.push_back(result_0_2);
      result.push_back(result_0_3);

      result.push_back(result_1_0);
      result.push_back(result_1_1);
      result.push_back(result_1_2);
      result.push_back(result_1_3);

      result.push_back(result_2_0);
      result.push_back(result_2_1);
      result.push_back(result_2_2);
      result.push_back(result_2_3);

      result.push_back(result_3_0);
      result.push_back(result_3_1);
      result.push_back(result_3_2);
      result.push_back(result_3_3);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return u->val[0] + v->val[0];
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsMatrixFormSemiImplicitOutlet : public WeakForm::MultiComponentMatrixFormSurf
  {
  public:
    EulerEquationsMatrixFormSemiImplicitOutlet(Hermes::vector<std::pair<unsigned int, 
                                                unsigned int> >coordinates, 
                                  std::string marker, double kappa) 
    : WeakForm::MultiComponentMatrixFormSurf(coordinates, marker), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { }

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext, 
               Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) {
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        double e_1_1[4] = {1, 0, 0, 0};
        double e_2_2[4] = {0, 1, 0, 0};
        double e_3_3[4] = {0, 0, 1, 0};
        double e_4_4[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];
      
        num_flux->P_plus(P_plus_1, w_L, e_1_1, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_2, w_L, e_2_2, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_3, w_L, e_3_3, e->nx[i], e->ny[i]);
        num_flux->P_plus(P_plus_4, w_L, e_4_4, e->nx[i], e->ny[i]);
        
        double e_1_1_2[4] = {1, 0, 0, 0};
        double e_2_2_2[4] = {0, 1, 0, 0};
        double e_3_3_2[4] = {0, 0, 1, 0};
        double e_4_4_2[4] = {0, 0, 0, 1};
              
        double P_minus_1[4];
        double P_minus_2[4];
        double P_minus_3[4];
        double P_minus_4[4];
      
        num_flux->P_minus(P_minus_1, w_L, e_1_1_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_2, w_L, e_2_2_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_3, w_L, e_3_3_2, e->nx[i], e->ny[i]);
        num_flux->P_minus(P_minus_4, w_L, e_4_4_2, e->nx[i], e->ny[i]);

        result_0_0 += wt[i] * (P_plus_1[0] + P_minus_1[0]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_1 += wt[i] * (P_plus_2[0] + P_minus_2[0]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_2 += wt[i] * (P_plus_3[0] + P_minus_3[0]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_3 += wt[i] * (P_plus_4[0] + P_minus_4[0]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_1_0 += wt[i] * (P_plus_1[1] + P_minus_1[1]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_1 += wt[i] * (P_plus_2[1] + P_minus_2[1]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_2 += wt[i] * (P_plus_3[1] + P_minus_3[1]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_3 += wt[i] * (P_plus_4[1] + P_minus_4[1]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_2_0 += wt[i] * (P_plus_1[2] + P_minus_1[2]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_1 += wt[i] * (P_plus_2[2] + P_minus_2[2]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_2 += wt[i] * (P_plus_3[2] + P_minus_3[2]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_3 += wt[i] * (P_plus_4[2] + P_minus_4[2]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_3_0 += wt[i] * (P_plus_1[3] + P_minus_1[3]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_1 += wt[i] * (P_plus_2[3] + P_minus_2[3]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_2 += wt[i] * (P_plus_3[3] + P_minus_3[3]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_3 += wt[i] * (P_plus_4[3] + P_minus_4[3]) * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
      }

      result.push_back(result_0_0);
      result.push_back(result_0_1);
      result.push_back(result_0_2);
      result.push_back(result_0_3);

      result.push_back(result_1_0);
      result.push_back(result_1_1);
      result.push_back(result_1_2);
      result.push_back(result_1_3);

      result.push_back(result_2_0);
      result.push_back(result_2_1);
      result.push_back(result_2_2);
      result.push_back(result_2_3);

      result.push_back(result_3_0);
      result.push_back(result_3_1);
      result.push_back(result_3_2);
      result.push_back(result_3_3);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return u->val[0] + v->val[0];
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsVectorFormSemiImplicitInlet : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    EulerEquationsVectorFormSemiImplicitInlet(Hermes::vector<unsigned int> coordinates, 
                                  std::string marker, double kappa) 
    : WeakForm::MultiComponentVectorFormSurf(coordinates, marker), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { }

    void value(int n, double *wt, Func<scalar> *u_ext[],
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext, 
               Hermes::vector<double>& result) const {
      double result_0 = 0;
      double result_1 = 0;
      double result_2 = 0;
      double result_3 = 0;

      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) {
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->rho_ext;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->v1_ext;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->rho_ext 
                 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->v2_ext;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->energy_ext;

        double P_minus[4];

        double w_B_temp[4];
        w_B_temp[0] = w_B[0];
        w_B_temp[1] = w_B[1];
        w_B_temp[2] = w_B[2];
        w_B_temp[3] = w_B[3];
      
        num_flux->P_minus(P_minus, w_B, w_B_temp, e->nx[i], e->ny[i]);

        result_0 += wt[i] * (P_minus[0]) * v->val[i];
        result_1 += wt[i] * (P_minus[1]) * v->val[i];
        result_2 += wt[i] * (P_minus[2]) * v->val[i];
        result_3 += wt[i] * (P_minus[3]) * v->val[i];
      }

      result.push_back(-result_0 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(-result_1 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(-result_2 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
      result.push_back(-result_3 * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau());
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0] + v->val[0];
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormTime : public WeakForm::MultiComponentVectorFormVol
  {
  public:
    EulerEquationsLinearFormTime(Hermes::vector<unsigned int> coordinates) 
         : WeakForm::MultiComponentVectorFormVol(coordinates) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[0], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[1], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[2], v));
      result.push_back(int_u_v<double, double>(n, wt, ext->fn[3], v));
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return int_u_v<Ord, Ord>(n, wt, v, v);
    }
  };

  class EulerEquationsMatrixFormSolidWall : public WeakForm::MultiComponentMatrixFormSurf
  {
  public:
    EulerEquationsMatrixFormSolidWall(Hermes::vector<std::pair<unsigned int, 
                                                unsigned int> >coordinates, 
                                      std::string marker, double kappa) 
         : WeakForm::MultiComponentMatrixFormSurf(coordinates, marker), kappa(kappa) {}

    void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, 
               ExtData<scalar> *ext, Hermes::vector<double>& result) const {
      double result_0_0 = 0;
      double result_0_1 = 0;
      double result_0_2 = 0;
      double result_0_3 = 0;

      double result_1_0 = 0;
      double result_1_1 = 0;
      double result_1_2 = 0;
      double result_1_3 = 0;

      double result_2_0 = 0;
      double result_2_1 = 0;
      double result_2_2 = 0;
      double result_2_3 = 0;

      double result_3_0 = 0;
      double result_3_1 = 0;
      double result_3_2 = 0;
      double result_3_3 = 0;

      for (int i = 0;i < n;i++) {
        double rho = ext->fn[0]->val[i];
        double v_1 = ext->fn[1]->val[i] / rho;
        double v_2 = ext->fn[2]->val[i] / rho;

        double P[4][4];
        for(unsigned int P_i = 0; P_i < 4; P_i++)
          for(unsigned int P_j = 0; P_j < 4; P_j++)
            P[P_i][P_j] = 0.0;


        P[1][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->nx[i] / 2;
        P[1][1] = (kappa - 1) * (-v_1) * e->nx[i];
        P[1][2] = (kappa - 1) * (-v_2) * e->nx[i];
        P[1][3] = (kappa - 1) * e->nx[i];

        P[2][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->ny[i] / 2;
        P[2][1] = (kappa - 1) * (-v_1) * e->ny[i];
        P[2][2] = (kappa - 1) * (-v_2) * e->ny[i];
        P[2][3] = (kappa - 1) * e->ny[i];

        result_0_0 += wt[i] * P[0][0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_1 += wt[i] * P[0][1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_2 += wt[i] * P[0][2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_0_3 += wt[i] * P[0][3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_1_0 += wt[i] * P[1][0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_1 += wt[i] * P[1][1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_2 += wt[i] * P[1][2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_1_3 += wt[i] * P[1][3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_2_0 += wt[i] * P[2][0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_1 += wt[i] * P[2][1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_2 += wt[i] * P[2][2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_2_3 += wt[i] * P[2][3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();

        result_3_0 += wt[i] * P[3][0] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_1 += wt[i] * P[3][1] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_2 += wt[i] * P[3][2] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
        result_3_3 += wt[i] * P[3][3] * u->val[i] * v->val[i] * static_cast<EulerEquationsWeakFormSemiImplicitMultiComponent*>(wf)->get_tau();
      }

      result.push_back(result_0_0);
      result.push_back(result_0_1);
      result.push_back(result_0_2);
      result.push_back(result_0_3);

      result.push_back(result_1_0);
      result.push_back(result_1_1);
      result.push_back(result_1_2);
      result.push_back(result_1_3);

      result.push_back(result_2_0);
      result.push_back(result_2_1);
      result.push_back(result_2_2);
      result.push_back(result_2_3);

      result.push_back(result_3_0);
      result.push_back(result_3_1);
      result.push_back(result_3_2);
      result.push_back(result_3_3);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return v->val[0] * Ord(2);
    }

    // Members.
    double kappa;
  };

  // Members.
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double pressure_ext;
  double energy_ext;
  double tau;
  EulerFluxes* euler_fluxes;
};

class EulerEquationsWeakFormExplicitCoupled : public EulerEquationsWeakFormExplicitMultiComponent
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicitCoupled(NumericalFlux* num_flux, double kappa,
                    double rho_ext, double v1_ext, double v2_ext,
                    double pressure_ext, std::string solid_wall_marker,
                    std::string inlet_marker, std::string outlet_marker,
                    Hermes::vector<std::string> natural_bc_concentration_markers,
                    Solution* prev_density, Solution* prev_density_vel_x,
                    Solution* prev_density_vel_y, Solution* prev_energy,
                    Solution* prev_concentration, double epsilon, bool fvm_only = false)
      : EulerEquationsWeakFormExplicitMultiComponent(num_flux, kappa, rho_ext, v1_ext, v2_ext, pressure_ext,
                                             solid_wall_marker, solid_wall_marker, inlet_marker,
                                             outlet_marker, prev_density, prev_density_vel_x,
                                             prev_density_vel_y, prev_energy, fvm_only, 5) {

    add_matrix_form(new EulerEquationsWeakFormExplicit::EulerEquationsBilinearFormTime(4));
    
    add_vector_form(new VectorFormConcentrationDiffusion(4, epsilon));
    vfvol.back()->ext.push_back(prev_density);
    vfvol.back()->ext.push_back(prev_density_vel_x);
    vfvol.back()->ext.push_back(prev_density_vel_y);
    vfvol.back()->ext.push_back(prev_energy);
    vfvol.back()->ext.push_back(prev_concentration);

    add_vector_form(new VectorFormConcentrationAdvection(4));
    vfvol.back()->ext.push_back(prev_density);
    vfvol.back()->ext.push_back(prev_density_vel_x);
    vfvol.back()->ext.push_back(prev_density_vel_y);
    vfvol.back()->ext.push_back(prev_energy);
    vfvol.back()->ext.push_back(prev_concentration);

    for(unsigned int i = 0;i < natural_bc_concentration_markers.size();i++) {
      add_vector_form_surf(new VectorFormConcentrationNatural(4, natural_bc_concentration_markers[i]));
      vfsurf.back()->ext.push_back(prev_density);
      vfsurf.back()->ext.push_back(prev_density_vel_x);
      vfsurf.back()->ext.push_back(prev_density_vel_y);
      vfsurf.back()->ext.push_back(prev_energy);
      vfsurf.back()->ext.push_back(prev_concentration);
    }
    
    add_vector_form_surf(new VectorFormConcentrationInterface(4));
    vfsurf.back()->ext.push_back(prev_density);
    vfsurf.back()->ext.push_back(prev_density_vel_x);
    vfsurf.back()->ext.push_back(prev_density_vel_y);
    vfsurf.back()->ext.push_back(prev_energy);
    vfsurf.back()->ext.push_back(prev_concentration);

    EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime(4);
    vector_form_time->ext.push_back(prev_density);
    vector_form_time->ext.push_back(prev_density_vel_x);
    vector_form_time->ext.push_back(prev_density_vel_y);
    vector_form_time->ext.push_back(prev_energy);
    vector_form_time->ext.push_back(prev_concentration);
    add_vector_form(vector_form_time);
  };

  // Destructor.
  ~EulerEquationsWeakFormExplicitCoupled() {};
protected:
  class VectorFormConcentrationDiffusion : public WeakForm::VectorFormVol
  {
  public:
    VectorFormConcentrationDiffusion(int i, double epsilon) 
          : WeakForm::VectorFormVol(i), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* concentration_prev = ext->fn[4];
      return - epsilon * int_grad_u_grad_v<Real, Scalar>(n, wt, concentration_prev, v) 
             * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }

    // Member.
    double epsilon;
  };

  class VectorFormConcentrationAdvection : public WeakForm::VectorFormVol
  {
  public:
    VectorFormConcentrationAdvection(int i) : WeakForm::VectorFormVol(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];
      Func<Real>* concentration_prev = ext->fn[4];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += wt[i] * concentration_prev->val[i] * ((density_vel_x_prev->val[i] 
                  * v->dx[i]) + (density_vel_y_prev->val[i] * v->dy[i]))
                  / density_prev->val[i];
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }
  };

  class VectorFormConcentrationNatural : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormConcentrationNatural(int i, std::string marker) 
          : WeakForm::VectorFormSurf(i, marker) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];
      Func<Real>* concentration_prev = ext->fn[4];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += wt[i] * v->val[i] * concentration_prev->val[i] 
                  * (density_vel_x_prev->val[i] * e->nx[i] + density_vel_y_prev->val[i] * e->ny[i])
                  / density_prev->val[i];
        // (OR: for inlet/outlet) result += wt[i] * v->val[i] * concentration_prev->val[i] 
        //      * (V1_EXT * e->nx[i] + V2_EXT * e->ny[i]);
      return - result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }
  };

  class VectorFormConcentrationInterface : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormConcentrationInterface(int i) : WeakForm::VectorFormSurf(i, H2D_DG_INNER_EDGE) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];
      Func<Real>* concentration_prev = ext->fn[4];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += 0.5 * wt[i] * (v->get_val_central(i) - v->get_val_neighbor(i)) *
        (
          (
            density_vel_x_prev->get_val_central(i) * concentration_prev->get_val_central(i)
            / density_prev->get_val_central(i)
            +
            density_vel_x_prev->get_val_neighbor(i) * concentration_prev->get_val_neighbor(i)
            / density_prev->get_val_neighbor(i)
          )
        );

      return - result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                 ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }
  };
};

class EulerEquationsWeakFormSemiImplicitCoupled : public EulerEquationsWeakFormSemiImplicitMultiComponent
{
public:
  // Constructors.
  EulerEquationsWeakFormSemiImplicitCoupled(NumericalFlux* num_flux, double kappa,
                    double rho_ext, double v1_ext, double v2_ext,
                    double pressure_ext, std::string solid_wall_marker,
                    std::string inlet_marker, std::string outlet_marker,
                    Hermes::vector<std::string> natural_bc_concentration_markers,
                    Solution* prev_density, Solution* prev_density_vel_x,
                    Solution* prev_density_vel_y, Solution* prev_energy,
                    Solution* prev_concentration, double epsilon, bool fvm_only = false)
      : EulerEquationsWeakFormSemiImplicitMultiComponent(num_flux, kappa, rho_ext, v1_ext, v2_ext, pressure_ext,
                                             solid_wall_marker, solid_wall_marker, inlet_marker,
                                             outlet_marker, prev_density, prev_density_vel_x,
                                             prev_density_vel_y, prev_energy, fvm_only, 5) {

    add_matrix_form(new EulerEquationsWeakFormExplicit::EulerEquationsBilinearFormTime(4));
    
    add_matrix_form(new MatrixFormConcentrationDiffusion(4, 4, epsilon));
    mfvol.back()->ext.push_back(prev_density);
    mfvol.back()->ext.push_back(prev_density_vel_x);
    mfvol.back()->ext.push_back(prev_density_vel_y);
    mfvol.back()->ext.push_back(prev_energy);

    add_matrix_form(new MatrixFormConcentrationAdvection(4, 4));
    mfvol.back()->ext.push_back(prev_density);
    mfvol.back()->ext.push_back(prev_density_vel_x);
    mfvol.back()->ext.push_back(prev_density_vel_y);
    mfvol.back()->ext.push_back(prev_energy);

    for(unsigned int i = 0;i < natural_bc_concentration_markers.size();i++) {
      add_matrix_form_surf(new MatrixFormConcentrationNatural(4, 4, natural_bc_concentration_markers[i]));
      mfsurf.back()->ext.push_back(prev_density);
      mfsurf.back()->ext.push_back(prev_density_vel_x);
      mfsurf.back()->ext.push_back(prev_density_vel_y);
      mfsurf.back()->ext.push_back(prev_energy);
    }

    EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime(4);
    vector_form_time->ext.push_back(prev_density);
    vector_form_time->ext.push_back(prev_density_vel_x);
    vector_form_time->ext.push_back(prev_density_vel_y);
    vector_form_time->ext.push_back(prev_energy);
    vector_form_time->ext.push_back(prev_concentration);
    add_vector_form(vector_form_time);
  };

  EulerEquationsWeakFormSemiImplicitCoupled(NumericalFlux* num_flux, double kappa,
                    double rho_ext, double v1_ext, double v2_ext,
                    double pressure_ext, std::string solid_wall_marker_bottom, std::string solid_wall_marker_top,
                    std::string inlet_marker, std::string outlet_marker,
                    Hermes::vector<std::string> natural_bc_concentration_markers,
                    Solution* prev_density, Solution* prev_density_vel_x,
                    Solution* prev_density_vel_y, Solution* prev_energy,
                    Solution* prev_concentration, double epsilon, bool fvm_only = false)
      : EulerEquationsWeakFormSemiImplicitMultiComponent(num_flux, kappa, rho_ext, v1_ext, v2_ext, pressure_ext,
                                             solid_wall_marker_bottom, solid_wall_marker_top, inlet_marker,
                                             outlet_marker, prev_density, prev_density_vel_x,
                                             prev_density_vel_y, prev_energy, fvm_only, 5) {

    add_matrix_form(new EulerEquationsWeakFormExplicit::EulerEquationsBilinearFormTime(4));
    
    add_matrix_form(new MatrixFormConcentrationAdvectionDiffusion(4, 4, epsilon));
    mfvol.back()->ext.push_back(prev_density);
    mfvol.back()->ext.push_back(prev_density_vel_x);
    mfvol.back()->ext.push_back(prev_density_vel_y);
    mfvol.back()->ext.push_back(prev_energy);

    /*
    add_matrix_form(new MatrixFormConcentrationAdvection(4, 4));
    mfvol.back()->ext.push_back(prev_density);
    mfvol.back()->ext.push_back(prev_density_vel_x);
    mfvol.back()->ext.push_back(prev_density_vel_y);
    mfvol.back()->ext.push_back(prev_energy);
    */

    for(unsigned int i = 0;i < natural_bc_concentration_markers.size();i++) {
      add_matrix_form_surf(new MatrixFormConcentrationNatural(4, 4, natural_bc_concentration_markers[i]));
      mfsurf.back()->ext.push_back(prev_density);
      mfsurf.back()->ext.push_back(prev_density_vel_x);
      mfsurf.back()->ext.push_back(prev_density_vel_y);
      mfsurf.back()->ext.push_back(prev_energy);
    }
    
    /*
    add_matrix_form_surf(new MatrixFormConcentrationInterface(4, 4));
    mfsurf.back()->ext.push_back(prev_density);
    mfsurf.back()->ext.push_back(prev_density_vel_x);
    mfsurf.back()->ext.push_back(prev_density_vel_y);
    mfsurf.back()->ext.push_back(prev_energy);
    */

    EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime(4);
    vector_form_time->ext.push_back(prev_density);
    vector_form_time->ext.push_back(prev_density_vel_x);
    vector_form_time->ext.push_back(prev_density_vel_y);
    vector_form_time->ext.push_back(prev_energy);
    vector_form_time->ext.push_back(prev_concentration);
    add_vector_form(vector_form_time);
  };

  // Destructor.
  ~EulerEquationsWeakFormSemiImplicitCoupled() {};
protected:

  class MatrixFormConcentrationAdvectionDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormConcentrationAdvectionDiffusion(int i, int j, double epsilon) 
          : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) {}
    
    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Real h_e = e->diam;
      Real s_c = 0.9;

      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      for (int i=0; i < n; i++) {
        Scalar v_1 = density_vel_x_prev->val[i] / density_prev->val[i];
        Scalar v_2 = density_vel_y_prev->val[i] / density_prev->val[i];


        result += wt[i] * (epsilon * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
                                   - (v_1 * u->val[i] * v->dx[i] + v_2 * u->val[i] * v->dy[i]));
      
        Real R_squared = pow(v_1 * u->dx[i] + v_2 * u->dy[i], 2.);
        Real R = sqrt(R_squared); //This just does fabs(b1 * u->dx[i] + b2 * u->dy[i]); but it can be parsed
        result += wt[i] * s_c * 0.5 * h_e * R * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) / (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
        
        Scalar b_norm = sqrt(v_1 * v_1 + v_2 * v_2);
        Real tau = 1. / sqrt( 9 * pow(4 * epsilon / pow(h_e, 2), 2) + pow(2 * b_norm / h_e, 2));
        result += wt[i] * tau * (-v_1 * v->dx[i] - v_2 * v->dy[i] + epsilon * v->laplace[i]) * (-v_1 * u->dx[i] - v_2 * u->dy[i] + epsilon * u->laplace[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double epsilon;
  };

  class MatrixFormConcentrationNatural : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormConcentrationNatural(int i, int j, std::string marker) 
          : WeakForm::MatrixFormSurf(i, j, marker) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += wt[i] * v->val[i] * u->val[i]
                  * (density_vel_x_prev->val[i] * e->nx[i] + density_vel_y_prev->val[i] * e->ny[i])
                  / density_prev->val[i];
        // (OR: for inlet/outlet) result += wt[i] * v->val[i] * concentration_prev->val[i] 
        //      * (V1_EXT * e->nx[i] + V2_EXT * e->ny[i]);
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();

    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
  
  class MatrixFormConcentrationAdvection : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormConcentrationAdvection(int i, int j) : WeakForm::MatrixFormVol(i, j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += wt[i] * (u->dx[i] * density_vel_x_prev->val[i] * v->val[i] / density_prev->val[i] + u->dy[i] * density_vel_y_prev->val[i] * v->val[i] / density_prev->val[i]);

      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext,u , v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }
  };

  class MatrixFormConcentrationDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormConcentrationDiffusion(int i, int j, double epsilon) 
          : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      return -epsilon * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) 
             * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double epsilon;
  };
  
  /*
  class MatrixFormConcentrationAdvection : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormConcentrationAdvection(int i, int j) : WeakForm::MatrixFormVol(i, j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      Scalar result = 0;
      for (int i = 0;i < n;i++)
        result += wt[i] * u->val[i] * ((density_vel_x_prev->val[i] * v->dx[i]) + (density_vel_y_prev->val[i] * v->dy[i])) / density_prev->val[i];

      return -result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext,u , v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return Ord(5);
    }
  };
  */
};
