#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "euler_util.h"

class EulerEquationsWeakFormStabilization : public WeakForm<double>
{
public:
  EulerEquationsWeakFormStabilization(MeshFunctionSharedPtr<double> prev_rho) : WeakForm<double>()
  {
    this->set_ext(prev_rho);
    add_vector_form_DG(new DGVectorFormIndicator());
  }

  class DGVectorFormIndicator : public VectorFormDG<double>
  {
  public:
    DGVectorFormIndicator() 
      : VectorFormDG<double>(0)
    {
    }

    double value(int n, double *wt, DiscontinuousFunc<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, DiscontinuousFunc<double>** ext) const 
    {
      double result = 0;
      
      for (int i = 0;i < n;i++)
        result += wt[i] * v->val[i] * (ext[0]->val[i] - ext[0]->val_neighbor[i]) * (ext[0]->val[i] - ext[0]->val_neighbor[i]);

      return result / (e->get_diam_approximation(n) * std::pow(e->get_area(n, wt), 0.75));
    }

    Ord ord(int n, double *wt, DiscontinuousFunc<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return v->val[0] * v->val[0] * Ord(6);
    }

    VectorFormDG<double>* clone() const { return new DGVectorFormIndicator; }
  };
};

class EulerEquationsWeakFormSemiImplicit : public WeakForm<double>
{
public:
  double kappa;
  bool fvm_only;
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> inlet_markers;
  Hermes::vector<std::string> outlet_markers;

  MeshFunctionSharedPtr<double> prev_density;
  MeshFunctionSharedPtr<double> prev_density_vel_x;
  MeshFunctionSharedPtr<double> prev_density_vel_y;
  MeshFunctionSharedPtr<double> prev_energy;

  // External state.
  Hermes::vector<double> rho_ext;
  Hermes::vector<double> v1_ext;
  Hermes::vector<double> v2_ext;
  Hermes::vector<double> pressure_ext;
  Hermes::vector<double> energy_ext;

  // Fluxes for calculation.
  EulerFluxes* euler_fluxes;

  // Discrete indicator in the case of Feistauer limiting.
  bool* discreteIndicator;
  int discreteIndicatorSize;

  // For cache handling.
  class EulerEquationsMatrixFormSurfSemiImplicit;
  class EulerEquationsMatrixFormSemiImplicitInletOutlet;
  bool cacheReadyDG;
  bool cacheReadySurf;
  double** P_plus_cache_DG;
  double** P_minus_cache_DG;
  double** P_plus_cache_surf;
  double** P_minus_cache_surf;

  // Utility.
  bool oneInflow;

  // Constructor for one inflow.
  EulerEquationsWeakFormSemiImplicit(double kappa, 
    double rho_ext, double v1_ext, double v2_ext, double pressure_ext,
    Hermes::vector<std::string> solid_wall_markers, Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers, 
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y,  MeshFunctionSharedPtr<double> prev_energy, 
    bool fvm_only = false, int num_of_equations = 4) :

  WeakForm<double>(num_of_equations), 
    kappa(kappa), 
    solid_wall_markers(solid_wall_markers), inlet_markers(inlet_markers), outlet_markers(outlet_markers), 
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), 
    fvm_only(fvm_only), 
    euler_fluxes(new EulerFluxes(kappa)), discreteIndicator(NULL)
  {
    oneInflow = true;

    this->rho_ext.push_back(rho_ext);
    this->v1_ext.push_back(v1_ext);
    this->v2_ext.push_back(v2_ext);
    this->pressure_ext.push_back(pressure_ext); 
    energy_ext.push_back(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa));

    P_plus_cache_DG = new double*[13];
    P_minus_cache_DG = new double*[13];
    P_plus_cache_surf = new double*[13];
    P_minus_cache_surf = new double*[13];

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
    {
      P_plus_cache_DG[coordinate_i] = new double[16];
      P_minus_cache_DG[coordinate_i] = new double[16];
      P_plus_cache_surf[coordinate_i] = new double[16];
      P_minus_cache_surf[coordinate_i] = new double[16];
    }

    for(int form_i = 0; form_i < 4; form_i++)
    {
      add_matrix_form(new EulerEquationsBilinearFormTime(form_i));

      add_vector_form(new EulerEquationsLinearFormTime(form_i));

      add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(form_i, rho_ext, v1_ext, v2_ext, energy_ext[0], inlet_markers, kappa));

      if(outlet_markers.size() > 0)
        add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(form_i, rho_ext, v1_ext, v2_ext, energy_ext[0], outlet_markers, kappa));

      for(int form_j = 0; form_j < 4; form_j++)
      {
        if(!fvm_only) 
          add_matrix_form(new EulerEquationsBilinearForm(form_i, form_j, euler_fluxes));

        EulerEquationsMatrixFormSurfSemiImplicit* formDG = new EulerEquationsMatrixFormSurfSemiImplicit(form_i, form_j, kappa, euler_fluxes, &this->cacheReadyDG, this->P_plus_cache_DG, this->P_minus_cache_DG);
        add_matrix_form_DG(formDG);

        EulerEquationsMatrixFormSemiImplicitInletOutlet* formSurf = new EulerEquationsMatrixFormSemiImplicitInletOutlet(form_i, form_j, rho_ext, v1_ext, v2_ext, energy_ext[0], inlet_markers, kappa, &this->cacheReadySurf, this->P_plus_cache_surf, this->P_minus_cache_surf);
        add_matrix_form_surf(formSurf);

        if(outlet_markers.size() > 0)
        {
          formSurf = new EulerEquationsMatrixFormSemiImplicitInletOutlet(form_i, form_j, rho_ext, v1_ext, v2_ext, energy_ext[0], outlet_markers, kappa, &this->cacheReadySurf, this->P_plus_cache_surf, this->P_minus_cache_surf);
          add_matrix_form_surf(formSurf);
        }

        add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(form_i, form_j, solid_wall_markers, kappa));
      }
    }

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  // Constructor for more inflows.
  EulerEquationsWeakFormSemiImplicit(double kappa, 
    Hermes::vector<double> rho_ext, Hermes::vector<double> v1_ext, Hermes::vector<double> v2_ext, Hermes::vector<double> pressure_ext,
    Hermes::vector<std::string> solid_wall_markers, Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers, 
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y,  MeshFunctionSharedPtr<double> prev_energy, 
    bool fvm_only = false, int num_of_equations = 4) :

  WeakForm<double>(num_of_equations), 
    kappa(kappa), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
    solid_wall_markers(solid_wall_markers), inlet_markers(inlet_markers), outlet_markers(outlet_markers), 
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), 
    fvm_only(fvm_only), 
    euler_fluxes(new EulerFluxes(kappa)), discreteIndicator(NULL)
  {
    oneInflow = false;
    
    for(unsigned int inlet_i = 0; inlet_i < inlet_markers.size(); inlet_i++)
      energy_ext.push_back(QuantityCalculator::calc_energy(rho_ext[inlet_i], rho_ext[inlet_i] * v1_ext[inlet_i], rho_ext[inlet_i] * v2_ext[inlet_i], pressure_ext[inlet_i], kappa));

    P_plus_cache_DG = new double*[13];
    P_minus_cache_DG = new double*[13];
    P_plus_cache_surf = new double*[13];
    P_minus_cache_surf = new double*[13];

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
    {
      P_plus_cache_DG[coordinate_i] = new double[16];
      P_minus_cache_DG[coordinate_i] = new double[16];
      P_plus_cache_surf[coordinate_i] = new double[16];
      P_minus_cache_surf[coordinate_i] = new double[16];
    }

    for(int form_i = 0; form_i < 4; form_i++)
    {
      add_matrix_form(new EulerEquationsBilinearFormTime(form_i));

      add_vector_form(new EulerEquationsLinearFormTime(form_i));

      for(unsigned int inlet_i = 0; inlet_i < inlet_markers.size(); inlet_i++)
        add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(form_i, rho_ext[inlet_i], v1_ext[inlet_i], v2_ext[inlet_i], energy_ext[inlet_i], inlet_markers[inlet_i], kappa));

      add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(form_i, 0, 0, 0, 0, outlet_markers, kappa));

      for(int form_j = 0; form_j < 4; form_j++)
      {
        if(!fvm_only) 
          add_matrix_form(new EulerEquationsBilinearForm(form_i, form_j, euler_fluxes));

        EulerEquationsMatrixFormSurfSemiImplicit* formDG = new EulerEquationsMatrixFormSurfSemiImplicit(form_i, form_j, kappa, euler_fluxes, &this->cacheReadyDG, this->P_plus_cache_DG, this->P_minus_cache_DG);
        add_matrix_form_DG(formDG);

        for(unsigned int inlet_i = 0; inlet_i < inlet_markers.size(); inlet_i++)
        {
          EulerEquationsMatrixFormSemiImplicitInletOutlet* formSurf = new EulerEquationsMatrixFormSemiImplicitInletOutlet(form_i, form_j, rho_ext[inlet_i], v1_ext[inlet_i], v2_ext[inlet_i], energy_ext[inlet_i], inlet_markers[inlet_i], kappa, &this->cacheReadySurf, this->P_plus_cache_surf, this->P_minus_cache_surf);
          add_matrix_form_surf(formSurf);
        }

        EulerEquationsMatrixFormSemiImplicitInletOutlet* formSurf = new EulerEquationsMatrixFormSemiImplicitInletOutlet(form_i, form_j, 0,0,0,0, outlet_markers, kappa, &this->cacheReadySurf, this->P_plus_cache_surf, this->P_minus_cache_surf);
        add_matrix_form_surf(formSurf);

        add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(form_i, form_j, solid_wall_markers, kappa));
      }
    }

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  virtual ~EulerEquationsWeakFormSemiImplicit()
  {
    delete this->euler_fluxes;

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
    {
      delete [] P_plus_cache_DG[coordinate_i];
      delete [] P_minus_cache_DG[coordinate_i];
      delete [] P_plus_cache_surf[coordinate_i];
      delete [] P_minus_cache_surf[coordinate_i];
    }

    delete [] P_plus_cache_DG;
    delete [] P_minus_cache_DG;
    delete [] P_plus_cache_surf;
    delete [] P_minus_cache_surf;
  }

  void set_active_edge_state(Element** e, int isurf)
  {
    this->cacheReadySurf = false;
  }

  void set_active_DG_state(Element** e, int isurf)
  {
    this->cacheReadyDG = false;
  }

  WeakForm<double>* clone() const
  {
    EulerEquationsWeakFormSemiImplicit* wf;
    if(this->oneInflow)
      wf = new EulerEquationsWeakFormSemiImplicit(this->kappa, this->rho_ext[0], this->v1_ext[0], this->v2_ext[0], this->pressure_ext[0], 
    this->solid_wall_markers, this->inlet_markers, this->outlet_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->fvm_only, this->neq);
    else
      wf = new EulerEquationsWeakFormSemiImplicit(this->kappa, this->rho_ext, this->v1_ext, this->v2_ext, this->pressure_ext, 
    this->solid_wall_markers, this->inlet_markers, this->outlet_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->fvm_only, this->neq);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      MeshFunctionSharedPtr<double> ext = this->ext[i]->clone();

      if(dynamic_cast<Solution<double>*>(this->ext[i].get()))
      {
        if((dynamic_cast<Solution<double>*>(this->ext[i].get()))->get_type() == HERMES_SLN)
          dynamic_cast<Solution<double>*>(ext.get())->set_type(HERMES_SLN);
      }
      wf->ext.push_back(ext);
    }

    wf->set_current_time_step(this->get_current_time_step());

    /*
    EulerEquationsWeakFormSemiImplicit* wf = new EulerEquationsWeakFormSemiImplicit(this->kappa, this->rho_ext, this->v1_ext, this->v2_ext, this->pressure_ext, 
    this->solid_wall_markers, this->inlet_markers, this->outlet_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->fvm_only, this->neq);


    if(this->discreteIndicator != NULL)
    {
    bool* discreteIndicatorLocal = new bool[this->discreteIndicatorSize];
    memcpy(discreteIndicatorLocal, this->discreteIndicator, this->discreteIndicatorSize * sizeof(bool));
    wf->set_discreteIndicator(discreteIndicatorLocal, this->discreteIndicatorSize);
    }
    */

    return wf;
  }

  void cloneMembers(const WeakForm<double>* otherWf)
  {
  }

  void set_stabilization(MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y, MeshFunctionSharedPtr<double> prev_energy, double nu_1, double nu_2) 
  {
    int mfvol_size = this->mfvol.size();
    int mfsurf_size = this->mfsurf.size();   

    add_matrix_form(new EulerEquationsFormStabilizationVol(0, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(1, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(2, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(3, nu_1));

    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 3, nu_2));

    for(unsigned int matrix_form_i = mfvol_size;matrix_form_i < this->mfvol.size();matrix_form_i++) 
    {
      mfvol.at(matrix_form_i)->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }

    for(unsigned int matrix_form_i = mfsurf_size;matrix_form_i < this->mfsurf.size();matrix_form_i++) 
    {
      mfsurf.at(matrix_form_i)->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }
  }

  void set_discreteIndicator(bool* discreteIndicator, int size)
  {
    this->discreteIndicator = discreteIndicator;
    this->discreteIndicatorSize = size;
  }

  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>* *ext) const 
    {
      return int_u_v<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return int_u_v<Ord, Ord>(n, wt, u, v);
    }

    MatrixFormVol<double>* clone() const { return new EulerEquationsBilinearFormTime(this->i); }
  };

  class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int i, int j, EulerFluxes* fluxes)
      : MatrixFormVol<double>(i, j), fluxes(fluxes) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, 
      Func<double>* *ext) const
    {
      double result = 0.;
      for (int point_i = 0; point_i < n;point_i++) 
      {
        double rho = ext[0]->val[point_i];
        double rho_v_x = ext[1]->val[point_i];
        double rho_v_y = ext[2]->val[point_i];
        double rho_e = ext[3]->val[point_i];

        switch(i)
        {
        case 0:
          switch(j)
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
          switch(j)
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
          switch(j)
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
          switch(j)
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

      return - result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    MatrixFormVol<double>* clone() const
    {
      return new EulerEquationsBilinearForm(this->i, this->j, this->fluxes);
    }

    EulerFluxes* fluxes;
  };

  class EulerEquationsMatrixFormSurfSemiImplicit : public MatrixFormDG<double>
  {
  public:
    EulerEquationsMatrixFormSurfSemiImplicit(int i, int j, double kappa, EulerFluxes* fluxes, bool* cacheReady, double** P_plus_cache, double** P_minus_cache) 
      : MatrixFormDG<double>(i, j), num_flux(new StegerWarmingNumericalFlux(kappa)), cacheReady(cacheReady), P_plus_cache(P_plus_cache), P_minus_cache(P_minus_cache), fluxes(fluxes) 
    {
    }

    ~EulerEquationsMatrixFormSurfSemiImplicit() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, 
      DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double w[4];
      double result = 0.;

      if(!(*this->cacheReady))
      {
        for (int point_i = 0; point_i < n; point_i++) 
        {
          {
            w[0] = ext[0]->val[point_i];
            w[1] = ext[1]->val[point_i];
            w[2] = ext[2]->val[point_i];
            w[3] = ext[3]->val[point_i];

            double e_1_1[4] = {1, 0, 0, 0};
            double e_2_1[4] = {0, 1, 0, 0};
            double e_3_1[4] = {0, 0, 1, 0};
            double e_4_1[4] = {0, 0, 0, 1};

            num_flux->P_plus(this->P_plus_cache[point_i], w, e_1_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 4, w, e_2_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 8, w, e_3_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 12, w, e_4_1, e->nx[point_i], e->ny[point_i]);

            w[0] = ext[0]->val_neighbor[point_i];
            w[1] = ext[1]->val_neighbor[point_i];
            w[2] = ext[2]->val_neighbor[point_i];
            w[3] = ext[3]->val_neighbor[point_i];

            double e_1_2[4] = {1, 0, 0, 0};
            double e_2_2[4] = {0, 1, 0, 0};
            double e_3_2[4] = {0, 0, 1, 0};
            double e_4_2[4] = {0, 0, 0, 1};

            num_flux->P_minus(this->P_minus_cache[point_i], w, e_1_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 4, w, e_2_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 8, w, e_3_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 12, w, e_4_2, e->nx[point_i], e->ny[point_i]);
          }
        }
        *(const_cast<EulerEquationsMatrixFormSurfSemiImplicit*>(this))->cacheReady = true;
      }

      int index = j * 4 + i;

      if(u->val == NULL)
        if(v->val == NULL)
          for (int point_i = 0; point_i < n; point_i++) 
            result -= wt[point_i] * (this->P_minus_cache[point_i][index] * u->val_neighbor[point_i]) * v->val_neighbor[point_i];
        else
          for (int point_i = 0; point_i < n; point_i++) 
            result += wt[point_i] * (this->P_minus_cache[point_i][index] * u->val_neighbor[point_i]) * v->val[point_i];
      else
        if(v->val == NULL)
          for (int point_i = 0; point_i < n; point_i++) 
            result -= wt[point_i] * (this->P_plus_cache[point_i][index] * u->val[point_i]) * v->val_neighbor[point_i];
        else
          for (int point_i = 0; point_i < n; point_i++) 
            result += wt[point_i] * (this->P_plus_cache[point_i][index] * u->val[point_i]) * v->val[point_i];
        
      return result * wf->get_current_time_step();
    }

    MatrixFormDG<double>* clone()  const
    { 
      EulerEquationsMatrixFormSurfSemiImplicit* form = new EulerEquationsMatrixFormSurfSemiImplicit(this->i, this->j, this->num_flux->kappa, this->fluxes, this->cacheReady, this->P_plus_cache, this->P_minus_cache);
      form->wf = this->wf;
      return form;
    }

    bool* cacheReady;
    double** P_plus_cache;
    double** P_minus_cache;
    StegerWarmingNumericalFlux* num_flux;
    EulerFluxes* fluxes;
  };

  class EulerEquationsMatrixFormSemiImplicitInletOutlet : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSemiImplicitInletOutlet(int i, int j, double rho_ext, double v1_ext, double v2_ext, double energy_ext, std::string marker, double kappa, bool* cacheReady, double** P_plus_cache, double** P_minus_cache) 
      : MatrixFormSurf<double>(i, j), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), energy_ext(energy_ext), num_flux(new StegerWarmingNumericalFlux(kappa)), cacheReady(cacheReady), P_plus_cache(P_plus_cache), P_minus_cache(P_minus_cache)
    { 
      set_area(marker);
    }
    EulerEquationsMatrixFormSemiImplicitInletOutlet(int i, int j, double rho_ext, double v1_ext, double v2_ext, double energy_ext, Hermes::vector<std::string> markers, double kappa, bool* cacheReady, double** P_plus_cache, double** P_minus_cache) 
      : MatrixFormSurf<double>(i, j), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), energy_ext(energy_ext), num_flux(new StegerWarmingNumericalFlux(kappa)), cacheReady(cacheReady), P_plus_cache(P_plus_cache), P_minus_cache(P_minus_cache)
    { 
      set_areas(markers); 
    }

    ~EulerEquationsMatrixFormSemiImplicitInletOutlet() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;

      if(!(*this->cacheReady))
      {
        for (int point_i = 0; point_i < n; point_i++) 
        {
          double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

          // Inner state.
          w_L[0] = ext[0]->val[point_i];
          w_L[1] = ext[1]->val[point_i];
          w_L[2] = ext[2]->val[point_i];
          w_L[3] = ext[3]->val[point_i];

          // Transformation of the inner state to the local coordinates.
          num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

          // Initialize the matrices.
          double T[4][4];
          double T_inv[4][4];
          for(unsigned int ai = 0; ai < 4; ai++) 
          {
            for(unsigned int aj = 0; aj < 4; aj++) 
            {
              T[ai][aj] = 0.0;
              T_inv[ai][aj] = 0.0;
            }
            alpha[ai] = 0;
            beta[ai] = 0;
            q_ji[ai] = 0;
            w_ji[ai] = 0;
            eigenvalues[ai] = 0;
          }

          // Calculate Lambda^-.
          num_flux->Lambda(eigenvalues);
          num_flux->T_1(T);
          num_flux->T_2(T);
          num_flux->T_3(T);
          num_flux->T_4(T);
          num_flux->T_inv_1(T_inv);
          num_flux->T_inv_2(T_inv);
          num_flux->T_inv_3(T_inv);
          num_flux->T_inv_4(T_inv);

          // "Prescribed" boundary state.
          w_B[0] = this->rho_ext;
          w_B[1] = this->rho_ext * this->v1_ext;
          w_B[2] = this->rho_ext * this->v2_ext;
          w_B[3] = this->energy_ext;

          num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

          for(unsigned int ai = 0; ai < 4; ai++)
            for(unsigned int aj = 0; aj < 4; aj++)
              alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

          for(unsigned int bi = 0; bi < 4; bi++)
            for(unsigned int bj = 0; bj < 4; bj++)
              beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

          for(unsigned int si = 0; si< 4; si++)
            for(unsigned int sj = 0; sj < 4; sj++)
              if(eigenvalues[sj] < 0)
                q_ji[si] += beta[sj] * T[si][sj];
              else
                q_ji[si] += alpha[sj] * T[si][sj];

          num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

          double w_temp[4];
          w_temp[0] = (w_ji[0] + w_L[0]) / 2;
          w_temp[1] = (w_ji[1] + w_L[1]) / 2;
          w_temp[2] = (w_ji[2] + w_L[2]) / 2;
          w_temp[3] = (w_ji[3] + w_L[3]) / 2;

          double e_1[4] = {1, 0, 0, 0};
          double e_2[4] = {0, 1, 0, 0};
          double e_3[4] = {0, 0, 1, 0};
          double e_4[4] = {0, 0, 0, 1};

          num_flux->P_plus(this->P_plus_cache[point_i], w_temp, e_1, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(this->P_plus_cache[point_i] + 4, w_temp, e_2, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(this->P_plus_cache[point_i] + 8, w_temp, e_3, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(this->P_plus_cache[point_i] + 12, w_temp, e_4, e->nx[point_i], e->ny[point_i]);
        }

        *(const_cast<EulerEquationsMatrixFormSemiImplicitInletOutlet*>(this))->cacheReady = true;
      }

      int index = j * 4 + i;
      for (int point_i = 0; point_i < n; point_i++) 
      {
        result += wt[point_i] * P_plus_cache[point_i][index] * u->val[point_i] * v->val[point_i];
      }

      return result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* clone()  const
    { 
      EulerEquationsMatrixFormSemiImplicitInletOutlet* form = new EulerEquationsMatrixFormSemiImplicitInletOutlet(this->i, this->j, this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->areas, this->num_flux->kappa, this->cacheReady, this->P_plus_cache, this->P_minus_cache);
      form->wf = this->wf;
      return form;
    }

    double rho_ext;
    double v1_ext;
    double v2_ext;
    double energy_ext;
    bool* cacheReady;
    double** P_plus_cache;
    double** P_minus_cache;
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsVectorFormSemiImplicitInletOutlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSemiImplicitInletOutlet(int i, double rho_ext, double v1_ext, double v2_ext, double energy_ext, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), energy_ext(energy_ext),
      num_flux(new StegerWarmingNumericalFlux(kappa)) 
    {
      set_area(marker);
    }
    EulerEquationsVectorFormSemiImplicitInletOutlet(int i, double rho_ext, double v1_ext, double v2_ext, double energy_ext, Hermes::vector<std::string> markers, double kappa) 
      : VectorFormSurf<double>(i), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), energy_ext(energy_ext),
      num_flux(new StegerWarmingNumericalFlux(kappa)) 
    {
      set_areas(markers);
    }

    ~EulerEquationsVectorFormSemiImplicitInletOutlet() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = this->rho_ext;
        w_B[1] = this->rho_ext * this->v1_ext;
        w_B[2] = this->rho_ext * this->v2_ext;
        w_B[3] = this->energy_ext;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        num_flux->P_minus(P_minus, w_temp, w_ji, e->nx[point_i], e->ny[point_i]);

        result += wt[point_i] * (P_minus[i]) * v->val[point_i];
      }

      return - result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* clone()  const
    { 
      EulerEquationsVectorFormSemiImplicitInletOutlet* form = new EulerEquationsVectorFormSemiImplicitInletOutlet(this->i, this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->areas, this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    double rho_ext;
    double v1_ext;
    double v2_ext;
    double energy_ext;
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const 
    {
      return int_u_v<double, double>(n, wt, ext[this->i], v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return int_u_v<Ord, Ord>(n, wt, ext[this->i], v);
    }

    VectorFormVol<double>* clone() const { return new EulerEquationsLinearFormTime(this->i); }
  };

  class EulerEquationsMatrixFormSolidWall : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSolidWall(int i, int j, Hermes::vector<std::string> markers, double kappa)
      : MatrixFormSurf<double>(i, j), kappa(kappa) {set_areas(markers);}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;

      for (int point_i = 0; point_i < n; point_i++) 
      {
        double rho = ext[0]->val[point_i];
        double v_1 = ext[1]->val[point_i] / rho;
        double v_2 = ext[2]->val[point_i] / rho;

        double P[4][4];
        for(unsigned int P_i = 0; P_i < 4; P_i++)
          for(unsigned int P_j = 0; P_j < 4; P_j++)
            P[P_i][P_j] = 0.0;

        P[1][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->nx[point_i] / 2;
        P[1][1] = (kappa - 1) * (-v_1) * e->nx[point_i];
        P[1][2] = (kappa - 1) * (-v_2) * e->nx[point_i];
        P[1][3] = (kappa - 1) * e->nx[point_i];

        P[2][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->ny[point_i] / 2;
        P[2][1] = (kappa - 1) * (-v_1) * e->ny[point_i];
        P[2][2] = (kappa - 1) * (-v_2) * e->ny[point_i];
        P[2][3] = (kappa - 1) * e->ny[point_i];

        result += wt[point_i] * P[i][j] * u->val[point_i] * v->val[point_i];
      }

      return result * wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* clone()  const
    {
      EulerEquationsMatrixFormSolidWall* form = new EulerEquationsMatrixFormSolidWall(this->i, this->j, this->areas, this->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double kappa;
  };

  class EulerEquationsFormStabilizationVol : public MatrixFormVol<double>
  {
  public:
    EulerEquationsFormStabilizationVol(int i, double nu_1) : MatrixFormVol<double>(i, i), nu_1(nu_1) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0.;
      if(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->id]) 
        result = int_grad_u_grad_v<double, double>(n, wt, u, v);
      return result * nu_1 * e->get_diam_approximation(n);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }
    MatrixFormVol<double>* clone() const
    {
      EulerEquationsFormStabilizationVol* form = new EulerEquationsFormStabilizationVol(this->i, nu_1);
      return form;
    }
  private:
    double nu_1;
  };

  class EulerEquationsFormStabilizationSurf : public MatrixFormDG<double>
  {
  public:
    EulerEquationsFormStabilizationSurf(int i, int j, double nu_2) 
      : MatrixFormDG<double>(i, j), nu_2(nu_2) {}

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, 
      DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double result = 0.;

      if(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->id] && static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->get_neighbor_id()])
      {
        if(u->val == NULL)
          if(v->val == NULL)
            for (int point_i = 0; point_i < n; point_i++) 
              result += wt[i] * (- u->val_neighbor[point_i]) * (- v->val_neighbor[point_i]);
          else
            for (int point_i = 0; point_i < n; point_i++) 
              result += wt[i] * (- u->val_neighbor[point_i]) * ( v->val[point_i]);
        else
          if(v->val == NULL)
            for (int point_i = 0; point_i < n; point_i++) 
              result += wt[i] * (u->val[point_i]) * (- v->val_neighbor[point_i]);
          else
            for (int point_i = 0; point_i < n; point_i++) 
              result += wt[i] * (u->val[point_i]) * ( v->val[point_i]);
      }

      return result * nu_2;
    }

    MatrixFormDG<double>* clone()  const
    {
      EulerEquationsFormStabilizationSurf* form = new EulerEquationsFormStabilizationSurf(this->i, this->j, nu_2);
      form->wf = this->wf;
      return form;
    }

    double nu_2;
  };

  
};

class EulerEquationsWeakFormSemiImplicitCoupledWithHeat : public EulerEquationsWeakFormSemiImplicit
{
public:
  MeshFunctionSharedPtr<double> prev_temp;
  double lambda, c_p;

  EulerEquationsWeakFormSemiImplicitCoupledWithHeat(double kappa, 
    double rho_ext, double v1_ext, double v2_ext, double pressure_ext,
    Hermes::vector<std::string> solid_wall_markers, Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers, 
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x, MeshFunctionSharedPtr<double> prev_density_vel_y, MeshFunctionSharedPtr<double> prev_energy, MeshFunctionSharedPtr<double> prev_temp, double lambda, double c_p): EulerEquationsWeakFormSemiImplicit(kappa, rho_ext, v1_ext, v2_ext, pressure_ext, solid_wall_markers, inlet_markers, outlet_markers, prev_density,
    prev_density_vel_x, prev_density_vel_y, prev_energy, false, 5), prev_temp(prev_temp), lambda(lambda), c_p(c_p)
  {
    add_matrix_form(new HeatBilinearFormTime(4, c_p, lambda));

    add_vector_form(new HeatLinearFormTime(4, c_p));

    this->ext.push_back(prev_temp);
  }

  WeakForm<double>* clone() const
  {
    EulerEquationsWeakFormSemiImplicitCoupledWithHeat* wf;
    wf = new EulerEquationsWeakFormSemiImplicitCoupledWithHeat(this->kappa, this->rho_ext[0], this->v1_ext[0], this->v2_ext[0], this->pressure_ext[0], 
      this->solid_wall_markers, this->inlet_markers, this->outlet_markers, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->prev_temp, this->lambda, this->c_p);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      MeshFunctionSharedPtr<double> ext = this->ext[i]->clone();

      if(dynamic_cast<Solution<double>*>(this->ext[i].get()))
      {
        if((dynamic_cast<Solution<double>*>(this->ext[i].get()))->get_type() == HERMES_SLN)
          dynamic_cast<Solution<double>*>(ext.get())->set_type(HERMES_SLN);
      }
      wf->ext.push_back(ext);
    }

    wf->set_current_time_step(this->get_current_time_step());

    return wf;
  }

  class HeatBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    HeatBilinearFormTime(int i, double c_p, double lambda) : MatrixFormVol<double>(i, i), c_p(c_p), lambda(lambda) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>* *ext) const 
    {
      double result = 0.;
      for (int point_i = 0; point_i < n; point_i++)
      {
        double rho = ext[0]->val[point_i];
        result += wt[point_i] * u->val[point_i] * v->val[point_i] * rho * c_p / this->wf->get_current_time_step();
        result += wt[point_i] * (u->dx[point_i] * v->dx[point_i] + u->dy[point_i] * v->dy[point_i]) * lambda;
        result += wt[point_i] * (ext[1]->val[point_i] * u->dx[point_i] + ext[2]->val[point_i] * u->dy[point_i]) * v->val[point_i] * c_p;
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      Ord result = Ord(0);
      for (int point_i = 0; point_i < n; point_i++)
      {
        Ord rho = ext[0]->val[point_i];
        result += wt[point_i] * (u->val[point_i] * v->val[point_i] * rho - u->dx[point_i] * v->dx[point_i] - u->dy[point_i] * v->dy[point_i] - rho * rho * (ext[1]->val[point_i] * u->dx[point_i] + ext[2]->val[point_i] * u->dy[point_i]) * v->val[point_i]);
      }
      return result;
    }

    MatrixFormVol<double>* clone() const { return new HeatBilinearFormTime(*this); }

    double c_p, lambda;
  };

  class HeatLinearFormTime : public VectorFormVol<double>
  {
  public:
    HeatLinearFormTime(int i, double c_p) 
      : VectorFormVol<double>(i), c_p(c_p) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const 
    {
      double result = 0.;
      for (int point_i = 0; point_i < n; point_i++)
      {
        double rho = ext[0]->val[point_i];
        result += wt[i] * v->val[point_i] * rho * ext[this->i]->val[point_i];
      }
      return result * c_p / this->wf->get_current_time_step();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const 
    {
      return wt[0] * v->val[0] * ext[0]->val[0] * ext[this->i]->val[0];
    }

    VectorFormVol<double>* clone() const { return new HeatLinearFormTime(*this); }

    double c_p;
  };
};