#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "util.h"

class MHDWeakForm : public WeakForm<double>
{
public:
  double kappa;

  Hermes::vector<std::string> inlet_markers;
  Hermes::vector<std::string> outlet_markers;

  MeshFunctionSharedPtr<double> prev_density;
  MeshFunctionSharedPtr<double> prev_density_vel_x;
  MeshFunctionSharedPtr<double> prev_density_vel_y;
  MeshFunctionSharedPtr<double> prev_energy;
  ///TODO Zbyle 4 predchozi reseni

  // External state.
  double external_state[4];
  double pressure_ext;
  ///TODO Zbyle slozky externiho stavoveho vektoru

  // Fluxes for calculation.
  EulerFluxes* euler_fluxes;

  // Numerical flux.
  NumericalFlux* numerical_flux;

    // Constructor for one inflow.
    ///TODO Doplnit do kontruktoru parametry pro dalsi slozky stavoveho vektoru
    MHDWeakForm(double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext,
    Hermes::vector<std::string> inlet_markers, Hermes::vector<std::string> outlet_markers,
    NumericalFlux* numerical_flux,
    MeshFunctionSharedPtr<double> prev_density, MeshFunctionSharedPtr<double> prev_density_vel_x,
    MeshFunctionSharedPtr<double> prev_density_vel_y, MeshFunctionSharedPtr<double> prev_energy);

    virtual ~MHDWeakForm();

  ///TODO V teto metode staci patricne upravit zavolani kontruktoru, clonovani pole 'ext' je jiz nezavisle na poctu prvku v poli.
    WeakForm<double>* clone() const;

    void cloneMembers(const WeakForm<double>* otherWf);

  // Ukazka - bilinearni forma pro integral u*v
  ///TODO Zde je dulezite dat pozor ze z nejakych duvodu celou rovnici vynasobime casovym krokem - cili tady neni
  /// jedna / tau, ale pak jsou ostatni integraly tim tau nasobene.
  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, Func<double>* *ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
      Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;
  };

  /// Ukazka bilinearni formy volumetricke - vyuziti fluxu.
  class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int i, int j, EulerFluxes* fluxes);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    EulerFluxes* fluxes;
  };

  class EulerEquationsLinearFormInterface : public VectorFormDG<double>
  {
  public:
    EulerEquationsLinearFormInterface(int i, NumericalFlux* num_flux);

    double value(int n, double *wt, Func<double> *u_ext[], DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormDG<double>* clone() const;

    NumericalFlux* num_flux;
  }; 

  class EulerEquationsLinearFormInlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormInlet(int i, Hermes::vector<std::string> areas);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;
  };

  class EulerEquationsLinearFormOutlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormOutlet(int i, Hermes::vector<std::string> areas);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], DiscontinuousFunc<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
      Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const;
  };
};
