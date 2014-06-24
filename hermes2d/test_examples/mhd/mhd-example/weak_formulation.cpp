#include "hermes2d.h"

#include "weak_formulation.h"


MHDWeakForm::MHDWeakForm(double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext,
  Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers,
  NumericalFlux* numerical_flux,
  MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x,
  MeshFunctionSharedPtr<double> prev_density_vel_y, MeshFunctionSharedPtr<double> prev_energy) :

  ///TODO Zde tech "8" znamena pocet rovnic vysledne slabe formulace - cili "8" bude nas pripad
  WeakForm<double>(8),
  kappa(kappa), inlet_markers(inlet_markers), outlet_markers(outlet_markers),
  numerical_flux(numerical_flux),
  ///TODO Doplnit sem predani parametru pro dalsi slozky stavoveho vektoru
  prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy),
  euler_fluxes(new EulerFluxes(kappa))
{
    this->external_state[0] = rho_ext;
    this->external_state[1] = rho_ext * v1_ext;
    this->external_state[2] = rho_ext * v2_ext;
    this->pressure_ext = pressure_ext;
    this->external_state[3] = rho_ext * QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa);

    ///TODO Doplnit zde registraci forem.
    for (int form_i = 0; form_i < 4; form_i++)
    {
      add_matrix_form(new EulerEquationsBilinearFormTime(form_i));

      for (int form_j = 0; form_j < 4; form_j++)
        add_matrix_form(new EulerEquationsBilinearForm(form_i, form_j, euler_fluxes));

      add_vector_form(new EulerEquationsLinearFormTime(form_i));

      add_vector_form_DG(new EulerEquationsLinearFormInterface(form_i, numerical_flux));

      add_vector_form_surf(new EulerEquationsLinearFormInlet(form_i, inlet_markers));

      add_vector_form_surf(new EulerEquationsLinearFormOutlet(form_i, outlet_markers));

      /// ...

      /// ...

      /// ...
    }

    ///TODO Doplnit zde registraci reseni z predchoziho casoveho kroku pro zbyle slozky stavoveho vektoru.
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  }

MHDWeakForm::~MHDWeakForm()
{
}

///TODO V teto metode staci patricne upravit zavolani kontruktoru, clonovani pole 'ext' je jiz nezavisle na poctu prvku v poli.
WeakForm<double>* MHDWeakForm::clone() const
{
  MHDWeakForm* wf = new MHDWeakForm(this->kappa, this->external_state[0], this->external_state[1] / this->external_state[0], this->external_state[2] / this->external_state[0], this->pressure_ext,
    this->inlet_markers, this->outlet_markers, this->numerical_flux, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy);

  wf->ext.clear();

  for (unsigned int i = 0; i < this->ext.size(); i++)
  {
    MeshFunctionSharedPtr<double> ext = this->ext[i]->clone();

    if (dynamic_cast<Solution<double>*>(this->ext[i].get()))
    {
      if ((dynamic_cast<Solution<double>*>(this->ext[i].get()))->get_type() == HERMES_SLN)
        dynamic_cast<Solution<double>*>(ext.get())->set_type(HERMES_SLN);
    }
    wf->ext.push_back(ext);
  }

  wf->set_current_time_step(this->get_current_time_step());

  return wf;
}

void MHDWeakForm::cloneMembers(const WeakForm<double>* otherWf)
{
}

MHDWeakForm::EulerEquationsBilinearFormTime::EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i)
{
}

double MHDWeakForm::EulerEquationsBilinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
  Geom<double> *e, Func<double>* *ext) const
{
  return int_u_v<double, double>(n, wt, u, v);
}

Ord MHDWeakForm::EulerEquationsBilinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
  Func<Ord>* *ext) const
{
  return int_u_v<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* MHDWeakForm::EulerEquationsBilinearFormTime::clone() const { return new EulerEquationsBilinearFormTime(this->i); }

MHDWeakForm::EulerEquationsBilinearForm::EulerEquationsBilinearForm(int i, int j, EulerFluxes* fluxes)
: MatrixFormVol<double>(i, j), fluxes(fluxes) {}

double MHDWeakForm::EulerEquationsBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e,
  Func<double>* *ext) const
{
  double result = 0.;
  for (int point_i = 0; point_i < n; point_i++)
  {
    double rho = ext[0]->val[point_i];
    double rho_v_x = ext[1]->val[point_i];
    double rho_v_y = ext[2]->val[point_i];
    double rho_e = ext[3]->val[point_i];

    switch (i)
    {
    case 0:
      switch (j)
      {
      case 0:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 1:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 2:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 3:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_0_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_0_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      }
      break;
    case 1:
      switch (j)
      {
      case 0:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 1:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 2:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 3:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_1_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_1_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      }
      break;
    case 2:
      switch (j)
      {
      case 0:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_0(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_0(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 1:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_1(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_1(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 2:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_2(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_2(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      case 3:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_2_3(rho, rho_v_x, rho_v_y, 0) * v->dx[point_i] + fluxes->A_2_2_3(rho, rho_v_x, rho_v_y, 0) * v->dy[point_i]);
        break;
      }
      break;
    case 3:
      switch (j)
      {
      case 0:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_0(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_0(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
        break;
      case 1:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_1(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_1(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
        break;
      case 2:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_2(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_2(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
        break;
      case 3:
        result += wt[point_i] * u->val[point_i] * (fluxes->A_1_3_3(rho, rho_v_x, rho_v_y, rho_e) * v->dx[point_i] + fluxes->A_2_3_3(rho, rho_v_x, rho_v_y, rho_e) * v->dy[point_i]);
        break;
      }
    }
  }

  return -result * wf->get_current_time_step();
}

Ord MHDWeakForm::EulerEquationsBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return Ord(10);
}

MatrixFormVol<double>* MHDWeakForm::EulerEquationsBilinearForm::clone() const
{
  return new EulerEquationsBilinearForm(this->i, this->j, this->fluxes);
}

MHDWeakForm::EulerEquationsLinearFormInterface::EulerEquationsLinearFormInterface(int i, NumericalFlux* num_flux)
: VectorFormDG<double>(i), num_flux(num_flux)
{
}

double MHDWeakForm::EulerEquationsLinearFormInterface::value(int n, double *wt, Func<double> *u_ext[], DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0;
  double w_L[4], w_R[4];
  for (int i = 0; i < n; i++)
  {
    w_L[0] = ext[0]->val[i];
    w_L[1] = ext[1]->val[i];
    w_L[2] = ext[2]->val[i];
    w_L[3] = ext[3]->val[i];

    w_R[0] = ext[0]->val_neighbor[i];
    w_R[1] = ext[1]->val_neighbor[i];
    w_R[2] = ext[2]->val_neighbor[i];
    w_R[3] = ext[3]->val_neighbor[i];

    result += wt[i] * v->val[i] * num_flux->numerical_flux(this->i, w_L, w_R, e->nx[i], e->ny[i]);
  }
  return result * wf->get_current_time_step();
}

Ord MHDWeakForm::EulerEquationsLinearFormInterface::ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0];
}

VectorFormDG<double>* MHDWeakForm::EulerEquationsLinearFormInterface::clone() const
{
  return new EulerEquationsLinearFormInterface(*this);
}

MHDWeakForm::EulerEquationsLinearFormInlet::EulerEquationsLinearFormInlet(int i, Hermes::vector<std::string> areas)
: VectorFormSurf<double>(i)
{
  this->set_areas(areas);
}

double MHDWeakForm::EulerEquationsLinearFormInlet::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  MHDWeakForm* wf_mhd = static_cast<MHDWeakForm*>(this->wf);
  double external_state_value = wf_mhd->external_state[this->i];

  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * external_state_value;

  return result * wf->get_current_time_step();
}

Ord MHDWeakForm::EulerEquationsLinearFormInlet::ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0];
}

VectorFormSurf<double>* MHDWeakForm::EulerEquationsLinearFormInlet::clone() const
{
  return new EulerEquationsLinearFormInlet(*this);
}

MHDWeakForm::EulerEquationsLinearFormOutlet::EulerEquationsLinearFormOutlet(int i, Hermes::vector<std::string> areas)
: VectorFormSurf<double>(i)
{
  this->set_areas(areas);
}

double MHDWeakForm::EulerEquationsLinearFormOutlet::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;

  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * ext[this->i]->val[i];

  return result * wf->get_current_time_step();
}

Ord MHDWeakForm::EulerEquationsLinearFormOutlet::ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0];
}

VectorFormSurf<double>* MHDWeakForm::EulerEquationsLinearFormOutlet::clone() const
{
  return new EulerEquationsLinearFormOutlet(*this);
}

MHDWeakForm::EulerEquationsLinearFormTime::EulerEquationsLinearFormTime(int i)
: VectorFormVol<double>(i) {}

double MHDWeakForm::EulerEquationsLinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
  Func<double>* *ext) const
{
  return int_u_v<double, double>(n, wt, ext[this->i], v);
}

Ord MHDWeakForm::EulerEquationsLinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
  Func<Ord>* *ext) const
{
  return int_u_v<Ord, Ord>(n, wt, ext[this->i], v);
}

VectorFormVol<double>* MHDWeakForm::EulerEquationsLinearFormTime::clone() const
{
  return new EulerEquationsLinearFormTime(this->i);
}
