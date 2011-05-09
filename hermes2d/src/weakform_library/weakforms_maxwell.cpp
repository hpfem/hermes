#include "../hermes2d.h"

namespace WeakFormsMaxwell {
  DefaultJacobianMagnetostatics::DefaultJacobianMagnetostatics(int i, int j, std::string area, scalar const_coeff,
    CubicSpline* c_spline,
    SymFlag sym,
    GeomType gt,
    int order_increase)
    : WeakForm::MatrixFormVol(i, j, area, sym), idx_j(j), const_coeff(const_coeff), 
    spline_coeff(c_spline), gt(gt),
    order_increase(order_increase) 
  { 
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }
  DefaultJacobianMagnetostatics::DefaultJacobianMagnetostatics(int i, int j, Hermes::vector<std::string> areas,
    scalar const_coeff,
    CubicSpline* c_spline,
    SymFlag sym,
    GeomType gt,
    int order_increase)
    : WeakForm::MatrixFormVol(i, j, areas, sym), idx_j(j), const_coeff(const_coeff), 
    spline_coeff(c_spline), gt(gt),
    order_increase(order_increase) 
  { 
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  scalar DefaultJacobianMagnetostatics::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar planar_part = 0;
      scalar axisym_part = 0;
      for (int i = 0; i < n; i++) {
        scalar B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
        if (std::abs(B_i) > 1e-12) {
          planar_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
            * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) {
            axisym_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i / e->y[i]
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * u_ext[idx_j]->val[i] * v->dy[i];
          }
          else if (gt == HERMES_AXISYM_Y) {
            axisym_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i / e->x[i]
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * u_ext[idx_j]->val[i] * v->dx[i];
          }
        }
        planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i)
          * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        if (gt == HERMES_AXISYM_X) {
          axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->y[i]
          * u->val[i] * v->dy[i];
        }
        else if (gt == HERMES_AXISYM_Y) {
          axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->x[i]
          * u->val[i] * v->dx[i];
        }
      }

      return planar_part + axisym_part;
  }

  Ord DefaultJacobianMagnetostatics::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord planar_part = 0;
      for (int i = 0; i < n; i++) {
        Ord B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
        planar_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i
          * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
          * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
        planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i)
          * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }

      // This increase is for the axisymmetric part. We are not letting the
      // Ord class do it since it would automatically choose the highest order
      // due to the nonpolynomial 1/r term.
      return planar_part * Ord(order_increase);
  }

  WeakForm::MatrixFormVol* DefaultJacobianMagnetostatics::clone() {
    return new DefaultJacobianMagnetostatics(*this);
  }


  DefaultResidualMagnetostatics::DefaultResidualMagnetostatics(int i, std::string area, scalar const_coeff,
    CubicSpline* c_spline,
    GeomType gt,
    int order_increase)
    : WeakForm::VectorFormVol(i, area), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), 
    gt(gt), order_increase(order_increase) 
  { 
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  DefaultResidualMagnetostatics::DefaultResidualMagnetostatics(int i, Hermes::vector<std::string> areas, scalar const_coeff, 
    CubicSpline* c_spline,
    GeomType gt, int order_increase)
    : WeakForm::VectorFormVol(i, areas), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt),
    order_increase(order_increase) 
  { 
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  scalar DefaultResidualMagnetostatics::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<scalar> *ext) const {
      scalar planar_part = 0;
      scalar axisym_part = 0;
      for (int i = 0; i < n; i++) {
        scalar B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
        planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) *
          (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        if (gt == HERMES_AXISYM_X) axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->y[i]
          * u_ext[idx_i]->val[i] * v->dy[i];
        else if (gt == HERMES_AXISYM_Y) axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->x[i]
          * u_ext[idx_i]->val[i] * v->dx[i];
      }
      return planar_part + axisym_part;
  }

  Ord DefaultResidualMagnetostatics::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord planar_part = 0;
      for (int i = 0; i < n; i++) {
        Ord B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
        planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) *
          (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
      }
      return planar_part * Ord(order_increase);
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::VectorFormVol* DefaultResidualMagnetostatics::clone() {
    return new DefaultResidualMagnetostatics(*this);
  }
};


/* Default volumetric matrix form \int_{area} u->val[i] * ((vel_x - e->y[i] * vel_ang) *
v->dx[i] + (vel_y + e->x[i] * vel_ang) * v->dy[i])
vel_x, vel_y... velocity components
vel_ang... velocity angle
*/

/*
class DefaultLinearMagnetostaticsVelocity : public WeakForm::MatrixFormVol
{
public:
DefaultLinearMagnetostaticsVelocity(int i, int j, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
: WeakForm::MatrixFormVol(i, j), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

DefaultLinearMagnetostaticsVelocity(int i, int j, std::string area, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
: WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
Geom<double> *e, ExtData<scalar> *ext) const {
scalar result = 0;
for (int i = 0; i < n; i++)
result += wt[i] * u->val[i] * ((vel_x - e->y[i] * vel_ang) * v->dx[i] +
(vel_y + e->x[i] * vel_ang) * v->dy[i]);

return -gamma * result;
}

Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
Geom<Ord> *e, ExtData<Ord> *ext) const {
Ord result = 0;
for (int i = 0; i < n; i++)
result += wt[i] * u->val[i] * (v->dx[i] + v->dy[i]);

return result;
}

// This is to make the form usable in rk_time_step().
WeakForm::MatrixFormVol* clone() {
return new DefaultLinearMagnetostaticsVelocity(*this);
}

private:
double gamma, vel_x, vel_y, vel_ang;
};
*/


/* Default volumetric vector form \int_{area} coeff
\nabla u_ext[idx_i] \cdot \nabla v d\bfx
coeff... constant parameter
*/

/*
class DefaultResidualLinearMagnetostatics : public WeakForm::VectorFormVol
{
public:
DefaultResidualLinearMagnetostatics(int i, scalar coeff, GeomType gt,
int order_increase)
: WeakForm::VectorFormVol(i), idx_i(i), coeff(coeff), gt(gt), order_increase(order_increase) { }
DefaultResidualLinearMagnetostatics(int i, std::string area, scalar coeff,
GeomType gt, int order_increase)
: WeakForm::VectorFormVol(i, area), idx_i(i), coeff(coeff), gt(gt), order_increase(order_increase) { }

scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
Geom<double> *e, ExtData<scalar> *ext) const {
scalar planar_part = int_grad_u_ext_grad_v<double, scalar>(n, wt, u_ext[idx_i], v);
scalar axisym_part = 0;
if (gt == HERMES_AXISYM_X)
axisym_part = int_u_dvdy_over_y<double, scalar>(n, wt, u_ext[idx_i], v, e);
else if (gt == HERMES_AXISYM_Y)
axisym_part = int_u_dvdx_over_x<double, scalar>(n, wt, u_ext[idx_i], v, e);

return coeff * (planar_part + axisym_part);
}

Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
Geom<Ord> *e, ExtData<Ord> *ext) const {
Ord planar_part = int_grad_u_ext_grad_v<Ord, Ord>(n, wt, u_ext[idx_i], v);
return planar_part * Ord(order_increase);
}

// This is to make the form usable in rk_time_step().
WeakForm::VectorFormVol* clone() {
return new DefaultResidualLinearMagnetostatics(*this);
}

private:
int idx_i;
scalar coeff;
GeomType gt;
int order_increase;
};
*/

/* Default volumetric vector form \int_{area} spline_coeff(u_ext[0])
\nabla u_ext[idx_i] \cdot \nabla v d\bfx
spline_coeff... non-constant parameter given by a cubic spline
*/



/* Default volumetric vector form \int_{area} rem/perm * (- sin(rem_ang / 180.0 * M_PI) *
v->dx[i + cos(rem_ang / 180.0 * M_PI) * v->dy[i])
rem... remanent induction
rem_ang... remanent induction angle
per... permeability
*/

/*
class DefaultLinearMagnetostaticsRemanence : public WeakForm::VectorFormVol
{
public:
DefaultLinearMagnetostaticsRemanence(int i, double perm, double rem, double rem_ang, GeomType gt)
: WeakForm::VectorFormVol(i), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

DefaultLinearMagnetostaticsRemanence(int i, std::string area, double perm, double rem, double rem_ang, GeomType gt)
: WeakForm::VectorFormVol(i, area), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
Geom<double> *e, ExtData<scalar> *ext) const {
scalar result = 0;
for (int i = 0; i < n; i++)
result += wt[i] * rem/perm * (- sin(rem_ang / 180.0 * M_PI) * v->dx[i]
+ cos(rem_ang / 180.0 * M_PI) * v->dy[i]);

return (gt == HERMES_PLANAR ? -1.0 : 1.0) * result;
}

Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
Geom<Ord> *e, ExtData<Ord> *ext) const {
Ord result = 0;
for (int i = 0; i < n; i++)
result += wt[i] * (v->dx[i] + v->dy[i]);

return result;
}

// This is to make the form usable in rk_time_step().
WeakForm::VectorFormVol* clone() {
return new DefaultLinearMagnetostaticsRemanence(*this);
}

private:
double perm, rem, rem_ang;
GeomType gt;
};
*/
