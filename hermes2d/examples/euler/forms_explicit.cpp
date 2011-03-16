#include "hermes2d.h"
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "../euler-util.cpp"

class EulerEquationsWeakFormExplicit : public WeakForm
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicit(double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
  std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
  Solution* prev_density, Solution* prev_density_vel_x, Solution* prev_density_vel_y, Solution* prev_energy, int num_of_equations = 4) :
  WeakForm(num_of_equations), num_flux(kappa), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), energy_ext(calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)) {
    add_matrix_form(new EulerEquationsBilinearFormTime(0));
    add_matrix_form(new EulerEquationsBilinearFormTime(1));
    add_matrix_form(new EulerEquationsBilinearFormTime(2));
    add_matrix_form(new EulerEquationsBilinearFormTime(3));

    add_vector_form(new EulerEquationsLinearFormDensity());
    add_vector_form(new EulerEquationsLinearFormDensityVelX());
    add_vector_form(new EulerEquationsLinearFormDensityVelY());
    add_vector_form(new EulerEquationsLinearFormDensityEnergy());

    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    add_vector_form_surf(new EulerEquationsLinearFormInterface(0));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(1));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(2));
    add_vector_form_surf(new EulerEquationsLinearFormInterface(3));

    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_bottom_marker));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_bottom_marker));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_bottom_marker));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_bottom_marker));

    if(solid_wall_bottom_marker != solid_wall_top_marker) {
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_top_marker));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_top_marker));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_top_marker));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_top_marker));
    }

    add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(0, inlet_marker));
    add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(1, inlet_marker));
    add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(2, inlet_marker));
    add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(3, inlet_marker));

    if(inlet_marker != outlet_marker) {
      add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(0, outlet_marker));
      add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(1, outlet_marker));
      add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(2, outlet_marker));
      add_vector_form_surf(new EulerEquationsLinearFormInletOutlet(3, outlet_marker));
    }

    for(unsigned int vector_form_i = 0; vector_form_i < this->vfvol.size(); vector_form_i++) {
      vfvol.at(vector_form_i)->ext.push_back(prev_density);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfvol.at(vector_form_i)->ext.push_back(prev_energy);
    }

    for(unsigned int vector_form_i = 0; vector_form_i < this->vfsurf.size(); vector_form_i++) {
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
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class EulerEquationsLinearFormDensity : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensity() : WeakForm::VectorFormVol(0) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * ext->fn[1]->val[i] * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class EulerEquationsLinearFormDensityVelX : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensityVelX() : WeakForm::VectorFormVol(1) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Scalar>
    Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return - ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (num_flux.kappa - 1.) * 
             ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho));
    }

    template<typename Scalar>
    Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return 2 * (rho_v_x / rho) - (num_flux.kappa - 1.) * (rho_v_x / rho);
    }

    template<typename Scalar>
    Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return - (num_flux.kappa - 1.) * (rho_v_y / rho);;
    }

    template<typename Scalar>
    Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return num_flux.kappa - 1.;
    }
  };

  class EulerEquationsLinearFormDensityVelY : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensityVelY() : WeakForm::VectorFormVol(2) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Scalar>
    Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return - rho_v_x * rho_v_y / (rho * rho);
    }

    template<typename Scalar>
    Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return rho_v_y / rho;
    }

    template<typename Scalar>
    Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return rho_v_x / rho;
    }

    template<typename Scalar>
    Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return 0;
    }
  };

  class EulerEquationsLinearFormDensityEnergy : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormDensityEnergy() : WeakForm::VectorFormVol(3) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], 0) * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Scalar>
    Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return - (rho_v_x * energy) / (rho * rho) - (rho_v_x / (rho * rho)) * (num_flux.kappa - 1.) * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_x / rho) * (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho));
    }

    template<typename Scalar>
    Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return (energy / rho) + (1 / rho) * (num_flux.kappa - 1.) * ( energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho))) - (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x) / (rho * rho));
    }

    template<typename Scalar>
    Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return - (num_flux.kappa - 1.) * (rho_v_x * rho_v_y) / (rho * rho);
    }

    template<typename Scalar>
    Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
      return rho_v_x / rho + (num_flux.kappa - 1.) * (rho_v_x / rho);
    }
  };

  class EulerEquationsLinearFormTime : public WeakForm::VectorFormVol
  {
  public:
    EulerEquationsLinearFormTime(int i) : WeakForm::VectorFormVol(i), component_i(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_u_v<Real, Scalar>(n, wt, ext->fn[component_i], v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    int component_i;
  };

  class EulerEquationsLinearFormInterface : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormInterface(int i) : WeakForm::VectorFormSurf(i, H2D_DG_INNER_EDGE), element(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Scalar w_l[4], w_r[4];
      for (int i = 0; i < n; i++) {
        w_l[0] = ext->fn[0]->get_val_central(i);
        w_r[0] = ext->fn[0]->get_val_neighbor(i);
    
        w_l[1] = ext->fn[1]->get_val_central(i);
        w_r[1] = ext->fn[1]->get_val_neighbor(i);

        w_l[2] = ext->fn[2]->get_val_central(i);
        w_r[2] = ext->fn[2]->get_val_neighbor(i);

        w_l[3] = ext->fn[3]->get_val_central(i);
        w_r[3] = ext->fn[3]->get_val_neighbor(i);

        result -= wt[i] * v->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->num_flux.numerical_flux_i(element,w_l,w_r,e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return v->val[0];
    }

    // Member.
    int element;
  };

  class EulerEquationsLinearFormSolidWall : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormSolidWall(int i, std::string marker) : WeakForm::VectorFormSurf(i, marker), element(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Scalar w01, w11, w21, w31;
      for (int i = 0; i < n; i++) 
      {
        w01 = ext->fn[0]->val[i];
    
        w11 = ext->fn[1]->val[i];

        w21 = ext->fn[2]->val[i];

        w31 = ext->fn[3]->val[i];

        double p_b = calc_pressure(w01, w11, w21, w31, static_cast<EulerEquationsWeakFormExplicit*>(wf)->num_flux.kappa);
    
        double flux[4];
        double alpha = atan2(e->ny[i], e->nx[i]);
        double mat_rot_inv[4][4];
        double flux_local[4];
        flux_local[0] = 0;
        flux_local[1] = p_b;
        flux_local[2] = 0;
        flux_local[3] = 0;
        static_cast<EulerEquationsWeakFormExplicit*>(wf)->num_flux.T_rot(mat_rot_inv, -alpha);
        static_cast<EulerEquationsWeakFormExplicit*>(wf)->num_flux.dot_vector(flux, mat_rot_inv, flux_local);

        result -= wt[i] * v->val[i] * flux[element];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return v->val[0];
    }

    // Member.
    int element;
  };

  class EulerEquationsLinearFormInletOutlet : public WeakForm::VectorFormSurf
  {
  public:
    EulerEquationsLinearFormInletOutlet(int i, std::string marker) : WeakForm::VectorFormSurf(i, marker), element(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Scalar w_l[4], w_r[4];

      for (int i = 0; i < n; i++) 
      {
        // Left (inner) state from the previous time level solution.
        w_l[0] = ext->fn[0]->val[i];
        w_l[1] = ext->fn[1]->val[i];
        w_l[2] = ext->fn[2]->val[i];
        w_l[3] = ext->fn[3]->val[i];
    
        w_r[0] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext;
        w_r[1] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v1_ext;
        w_r[2] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v2_ext;
        w_r[3] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->energy_ext;
    
        /*
        // The inlet part (left part of the computational domain).
        if(e->nx[i] < 0)
        {
          // Boundary state calculation.
          double rho_b = bc_density(e->y[i]);
          double velocity_x_b = bc_density_vel_x(e->y[i]) / bc_density(e->y[i]);
          double velocity_y_b = bc_density_vel_y(e->y[i]) / bc_density(e->y[i]);

          // Sound speed on the left (inner) side of the boundary.
          double sound_speed_l = calc_sound_speed(w_l[0], w_l[1], w_l[2], w_l[3]);

          // Intersection state calculation (marked with an underscore1 (_1)).
          double sound_speed_1 = sound_speed_l + (num_flux.R/num_flux.c_v) * (w_l[1]/w_l[0] - velocity_x_b);
          double rho_1 = std::pow(sound_speed_1*sound_speed_1*w_l[0]/(num_flux.kappa*calc_pressure(w_l[0], w_l[1], w_l[2], w_l[3])), num_flux.c_v/num_flux.R) * w_l[0];
          double velocity_x_1 = velocity_x_b;
          double velocity_y_1 = w_l[2] / w_l[0];

          // Boundary pressure calculated from the intersection state.
          double p_b = rho_1 * sound_speed_1 * sound_speed_1 / num_flux.kappa;
          // Calculation of the energy component of the intersection state.
          double energy_1 = calc_energy<double>(rho_1, velocity_x_1* rho_1, velocity_y_1 * rho_1, p_b);

          // Calculation of the state for inflow/outlow velocities above the local speed of sound.
          double sound_speed_l_star = num_flux.R/(num_flux.c_v * (2+num_flux.R/num_flux.c_v)) * w_l[1] / w_l[0] + 2 * sound_speed_l / (2+num_flux.R/num_flux.c_v);
          double rho_l_star = std::pow(sound_speed_l_star/sound_speed_l, 2*num_flux.c_v / num_flux.R) * w_l[0];
          double velocity_x_l_star = sound_speed_l_star;
          double velocity_y_l_star = w_l[2] / w_l[0];
          double p_l_star = rho_l_star * sound_speed_l_star * sound_speed_l_star / num_flux.kappa;
          double energy_l_star = calc_energy<double>(rho_l_star, velocity_x_l_star * rho_l_star, velocity_y_l_star * rho_l_star, p_l_star);

          // Inflow velocity below the local speed of sound (of the intersection state).
          if(velocity_x_b < sound_speed_1)
          {
            //Ondrej's code.
            double alpha = atan2(e->ny[i], e->nx[i]);
            double mat_rot_inv[4][4];
            double flux_local[4];
            double flux_global[4];
            double mat_rot[4][4];
            num_flux.T_rot(mat_rot, alpha);
            flux_global[0] = rho_1;
            flux_global[1] = velocity_x_1 * rho_1;
            flux_global[2] = velocity_y_1 * rho_1;
            flux_global[3] = energy_1;
            num_flux.dot_vector(flux_local, mat_rot, flux_global);
            flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            num_flux.T_rot(mat_rot_inv, -alpha);
            num_flux.dot_vector(flux, mat_rot_inv, flux);
          }
          // Inflow velocity above the local speed of sound (of the intersection state).
          else
          {
            //Ondrej's code.
            double alpha = atan2(e->ny[i], e->nx[i]);
            double mat_rot_inv[4][4];
            double flux_local[4];
            double flux_global[4];
            double mat_rot[4][4];
            num_flux.T_rot(mat_rot, alpha);
            flux_global[0] = rho_l_star;
            flux_global[1] = velocity_x_l_star * rho_l_star;
            flux_global[2] = velocity_y_l_star * rho_l_star;
            flux_global[3] = energy_l_star;
            num_flux.dot_vector(flux_local, mat_rot, flux_global);
            flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            num_flux.T_rot(mat_rot_inv, -alpha);
            num_flux.dot_vector(flux, mat_rot_inv, flux);
          }
        }
        // The outlet part (the right part of the boundary).
        else
        {
          // These calculations are the same as above.
          double p_b = bc_pressure(e->y[i]);
          double rho_b = w_l[0] * std::pow(p_b/calc_pressure(w_l[0], w_l[1], w_l[2], w_l[3]), (1/num_flux.kappa));
          double velocity_x_b = (w_l[1] / w_l[0]) + 2*(num_flux.c_v/num_flux.R)*(calc_sound_speed<double>(w_l[0], w_l[1], w_l[2], w_l[3]) - std::sqrt(num_flux.kappa * p_b / rho_b));
          double velocity_y_b = w_l[2] / w_l[0];
          double energy_b = calc_energy<double>(rho_b, velocity_x_b*rho_b, velocity_y_b*rho_b, p_b);

          double sound_speed_l_star = num_flux.R/(num_flux.c_v * (2+num_flux.R/num_flux.c_v)) * w_l[1] / w_l[0] + 2 * calc_sound_speed<double>(w_l[0], w_l[1], w_l[2], w_l[3]) / (2+num_flux.R/num_flux.c_v);
          double rho_l_star = std::pow(sound_speed_l_star/calc_sound_speed<double>(w_l[0], w_l[1], w_l[2], w_l[3]), 2*num_flux.c_v / num_flux.R) * w_l[0];
          double velocity_x_l_star = sound_speed_l_star;
          double velocity_y_l_star = w_l[2] / w_l[0];
          double p_l_star = rho_l_star * sound_speed_l_star * sound_speed_l_star / num_flux.kappa;
          double energy_l_star = calc_energy<double>(rho_l_star, velocity_x_l_star * rho_l_star, velocity_y_l_star * rho_l_star, p_l_star);

          double sound_speed_b = calc_sound_speed(rho_b, velocity_x_b*rho_b, velocity_y_b*rho_b, energy_b);
          //  Inflow velocity below the local speed of sound (of the intersection state).
          if(velocity_x_b < sound_speed_b)
          {
            //Ondrej's code.
            double alpha = atan2(e->ny[i], e->nx[i]);
            double mat_rot_inv[4][4];
            double flux_local[4];
            double flux_global[4];
            double mat_rot[4][4];
            num_flux.T_rot(mat_rot, alpha);
            flux_global[0] = rho_b;
            flux_global[1] = velocity_x_b * rho_b;
            flux_global[2] = velocity_y_b * rho_b;
            flux_global[3] = energy_b;
            num_flux.dot_vector(flux_local, mat_rot, flux_global);
            flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            num_flux.T_rot(mat_rot_inv, -alpha);
            num_flux.dot_vector(flux, mat_rot_inv, flux);
          }
          //  Outflow velocity above the local speed of sound (of the intersection state).
          else
          {
            //Ondrej's code.
            double alpha = atan2(e->ny[i], e->nx[i]);
            double mat_rot_inv[4][4];
            double flux_local[4];
            double flux_global[4];
            double mat_rot[4][4];
            num_flux.T_rot(mat_rot, alpha);
            flux_global[0] = rho_l_star;
            flux_global[1] = velocity_x_l_star * rho_l_star;
            flux_global[2] = velocity_y_l_star * rho_l_star;
            flux_global[3] = energy_l_star;
            num_flux.dot_vector(flux_local, mat_rot, flux_global);
            flux[0] = f_x(0, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[1] = f_x(1, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[2] = f_x(2, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            flux[3] = f_x(3, flux_local[0], flux_local[1], flux_local[2], flux_local[3]);
            num_flux.T_rot(mat_rot_inv, -alpha);
            num_flux.dot_vector(flux, mat_rot_inv, flux);
          }
        }
        */
    
        result -= wt[i] * v->val[i] * static_cast<EulerEquationsWeakFormExplicit*>(wf)->num_flux.numerical_flux_i(element,w_l,w_r,e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return v->val[0];
    }

    // Member.
    int element;
  };


  // Eulerian fluxes.
  //////////////////////////////////////////////////
  ////////First flux////////////////////////////////
  ////////////////////////////////////////////////// 
  template<typename Scalar>
  Scalar A_1_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_1_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 1;
  }

  template<typename Scalar>
  Scalar A_1_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_1_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (num_flux.kappa - 1.) * 
           ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 2 * (rho_v_x / rho) - (num_flux.kappa - 1.) * (rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (num_flux.kappa - 1.) * (rho_v_y / rho);;
  }

  template<typename Scalar>
  Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return num_flux.kappa - 1.;
  }

  template<typename Scalar>
  Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - rho_v_x * rho_v_y / (rho * rho);
  }

  template<typename Scalar>
  Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_y / rho;
  }

  template<typename Scalar>
  Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_x / rho;
  }

  template<typename Scalar>
  Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (rho_v_x * energy) / (rho * rho) - (rho_v_x / (rho * rho)) * (num_flux.kappa - 1.) * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_x / rho) * (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return (energy / rho) + (1 / rho) * (num_flux.kappa - 1.) * ( energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho))) - (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (num_flux.kappa - 1.) * (rho_v_x * rho_v_y) / (rho * rho);
  }

  template<typename Scalar>
  Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_x / rho + (num_flux.kappa - 1.) * (rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar f_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return A_1_0_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
      + A_1_0_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
      + A_1_0_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
      + A_1_0_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
  }

  template<typename Scalar>
  Scalar f_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return A_1_1_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
      + A_1_1_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
      + A_1_1_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
      + A_1_1_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
  }

  template<typename Scalar>
  Scalar f_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return A_1_2_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
      + A_1_2_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
      + A_1_2_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
      + A_1_2_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
  }

  template<typename Scalar>
  Scalar f_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return A_1_3_0<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho 
      + A_1_3_1<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_x 
      + A_1_3_2<Scalar>(rho, rho_v_x, rho_v_y, energy) * rho_v_y
      + A_1_3_3<Scalar>(rho, rho_v_x, rho_v_y, energy) * energy;
  }

  double f_x(int i, double w0, double w1, double w3, double w4) {
    if(i == 0)
      return f_1_0<double>(w0, w1, w3, w4);
    if(i == 1)
      return f_1_1<double>(w0, w1, w3, w4);
    if(i == 2)
      return f_1_2<double>(w0, w1, w3, w4);
    if(i == 3)
      return f_1_3<double>(w0, w1, w3, w4);
    return 0.0;
  }

  /////////////////////////////////////////////////
  ////////Second flux jacobian/////////////////////
  /////////////////////////////////////////////////
  template<typename Scalar>
  Scalar A_2_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_2_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_2_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 1;
  }

  template<typename Scalar>
  Scalar A_2_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_2_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - rho_v_x * rho_v_y / (rho * rho);
  }

  template<typename Scalar>
  Scalar A_2_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_y / rho;
  }

  template<typename Scalar>
  Scalar A_2_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_x / rho;
  }

  template<typename Scalar>
  Scalar A_2_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 0;
  }

  template<typename Scalar>
  Scalar A_2_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (num_flux.kappa - 1.) * (rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_2_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return 2 * (rho_v_y / rho) - (num_flux.kappa - 1.) * (rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_2_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return num_flux.kappa - 1.;
  }

  template<typename Scalar>
  Scalar A_2_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (num_flux.kappa - 1.) * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) * (num_flux.kappa - 1.) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return - (num_flux.kappa - 1.) * (rho_v_x * rho_v_y) / (rho * rho);
  }

  template<typename Scalar>
  Scalar A_2_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return (energy / rho) + (1 / rho) * (num_flux.kappa - 1.) * ( energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho))) - (num_flux.kappa - 1.) * ((rho_v_y * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return rho_v_y / rho + (num_flux.kappa - 1.) * (rho_v_y / rho);
  }

  // Members.
  NumericalFlux num_flux;
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double energy_ext;
  double tau;
};

// The parameter variant in the constructor has the following meaning:
// 1 - Dirichlet condition (concentration production) on the inlet.
// 2 - Dirichlet condition (concentration production) on the bottom.
// 3 - Dirichlet condition (concentration production) on the top.
class EulerEquationsWeakFormExplicitCoupled : public EulerEquationsWeakFormExplicit
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicitCoupled(int variant, double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
  std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
  Solution* prev_density, Solution* prev_density_vel_x, Solution* prev_density_vel_y, Solution* prev_energy, Solution* prev_concentration, double epsilon) :
  EulerEquationsWeakFormExplicit(kappa, rho_ext, v1_ext, v2_ext, pressure_ext, solid_wall_bottom_marker, solid_wall_top_marker, inlet_marker, outlet_marker, prev_density,
  prev_density_vel_x, prev_density_vel_y, prev_energy, 5) {

    add_matrix_form(new EulerEquationsBilinearFormTime(4));
    
    add_vector_form(new VectorFormConcentrationDiffusion(4, epsilon));
    vfvol.back()->ext.push_back(prev_concentration);
    vfvol.back()->ext.push_back(prev_density);
    vfvol.back()->ext.push_back(prev_density_vel_x);
    vfvol.back()->ext.push_back(prev_density_vel_y);
    
    add_vector_form(new VectorFormConcentrationAdvection(4));
    vfvol.back()->ext.push_back(prev_concentration);
    vfvol.back()->ext.push_back(prev_density);
    vfvol.back()->ext.push_back(prev_density_vel_x);
    vfvol.back()->ext.push_back(prev_density_vel_y);

    if(variant != 1) {
      add_vector_form_surf(new VectorFormConcentrationNatural(4, inlet_marker));
      vfsurf.back()->ext.push_back(prev_concentration);
      vfsurf.back()->ext.push_back(prev_density);
      vfsurf.back()->ext.push_back(prev_density_vel_x);
      vfsurf.back()->ext.push_back(prev_density_vel_y);
    }
    if(variant != 2) {
      add_vector_form_surf(new VectorFormConcentrationNatural(4, solid_wall_bottom_marker));
      vfsurf.back()->ext.push_back(prev_concentration);
      vfsurf.back()->ext.push_back(prev_density);
      vfsurf.back()->ext.push_back(prev_density_vel_x);
      vfsurf.back()->ext.push_back(prev_density_vel_y);
    }
    if(variant != 3) {
      add_vector_form_surf(new VectorFormConcentrationNatural(4, solid_wall_top_marker));
      vfsurf.back()->ext.push_back(prev_concentration);
      vfsurf.back()->ext.push_back(prev_density);
      vfsurf.back()->ext.push_back(prev_density_vel_x);
      vfsurf.back()->ext.push_back(prev_density_vel_y);
    }

    add_vector_form_surf(new VectorFormConcentrationNatural(4, outlet_marker));
    vfsurf.back()->ext.push_back(prev_concentration);
    vfsurf.back()->ext.push_back(prev_density);
    vfsurf.back()->ext.push_back(prev_density_vel_x);
    vfsurf.back()->ext.push_back(prev_density_vel_y);

    add_vector_form_surf(new VectorFormConcentrationInterface(4));
    vfsurf.back()->ext.push_back(prev_concentration);
    vfsurf.back()->ext.push_back(prev_density);
    vfsurf.back()->ext.push_back(prev_density_vel_x);
    vfsurf.back()->ext.push_back(prev_density_vel_y);

    EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsLinearFormTime(4);
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
    VectorFormConcentrationDiffusion(int i, double epsilon) : WeakForm::VectorFormVol(i), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Real>* concentration_prev = ext->fn[0];
      return - epsilon * int_grad_u_grad_v<Real, Scalar>(n, wt, concentration_prev, v) * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    double epsilon;
  };

  class VectorFormConcentrationAdvection : public WeakForm::VectorFormVol
  {
  public:
    VectorFormConcentrationAdvection(int i) : WeakForm::VectorFormVol(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Real>* concentration_prev = ext->fn[0];
      Func<Real>* density_prev = ext->fn[1];
      Func<Real>* density_vel_x_prev = ext->fn[2];
      Func<Real>* density_vel_y_prev = ext->fn[3];

      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * concentration_prev->val[i] * ((density_vel_x_prev->val[i] * v->dx[i]) + (density_vel_y_prev->val[i] * v->dy[i]))
                  / density_prev->val[i];
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class VectorFormConcentrationNatural : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormConcentrationNatural(int i, std::string marker) : WeakForm::VectorFormSurf(i, marker) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Real>* concentration_prev = ext->fn[0];
      Func<Real>* density_prev = ext->fn[1];
      Func<Real>* density_vel_x_prev = ext->fn[2];
      Func<Real>* density_vel_y_prev = ext->fn[3];

      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * v->val[i] * concentration_prev->val[i] * (density_vel_x_prev->val[i] * e->nx[i] + density_vel_y_prev->val[i] * e->ny[i])
                  / density_prev->val[i];
        // (OR: for inlet/outlet) result += wt[i] * v->val[i] * concentration_prev->val[i] * (V1_EXT * e->nx[i] + V2_EXT * e->ny[i]);
      return - result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class VectorFormConcentrationInterface : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormConcentrationInterface(int i) : WeakForm::VectorFormSurf(i, H2D_DG_INNER_EDGE) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Real>* concentration_prev = ext->fn[0];
      Func<Real>* density_prev = ext->fn[1];
      Func<Real>* density_vel_x_prev = ext->fn[2];
      Func<Real>* density_vel_y_prev = ext->fn[3];

      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * v->val[i] *
                  (
                    (
                      density_vel_x_prev->get_val_central(i) * concentration_prev->get_val_central(i) / density_prev->get_val_central(i)
                      -
                      density_vel_x_prev->get_val_neighbor(i) * concentration_prev->get_val_neighbor(i) / density_prev->get_val_neighbor(i)
                    ) * e->nx[i]
                    + 
                    (
                      density_vel_y_prev->get_val_central(i) * concentration_prev->get_val_central(i) / density_prev->get_val_central(i)
                      -
                      density_vel_y_prev->get_val_neighbor(i) * concentration_prev->get_val_neighbor(i) / density_prev->get_val_neighbor(i)
                    ) * e->ny[i]
                  );
      return - result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(20);
    }
  };
};