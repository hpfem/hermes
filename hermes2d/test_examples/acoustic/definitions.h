#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

#pragma region forms

template<typename Scalar>
class volume_matrix_acoustic_transient_planar_linear_form_1_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_acoustic_transient_planar_linear_form_1_1(unsigned int i, unsigned int j, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  MatrixFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class volume_matrix_acoustic_transient_planar_linear_form_1_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_acoustic_transient_planar_linear_form_1_2(unsigned int i, unsigned int j, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  MatrixFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class volume_matrix_acoustic_transient_planar_linear_form_2_2 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_acoustic_transient_planar_linear_form_2_2(unsigned int i, unsigned int j, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  MatrixFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class volume_matrix_acoustic_transient_planar_linear_form_2_1 : public MatrixFormVol<Scalar>
{
public:
  volume_matrix_acoustic_transient_planar_linear_form_2_1(unsigned int i, unsigned int j, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
    Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
    Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  MatrixFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class volume_vector_acoustic_transient_planar_linear_form_1_2 : public VectorFormVol<Scalar>
{
public:
  volume_vector_acoustic_transient_planar_linear_form_1_2(unsigned int i, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  VectorFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class volume_vector_acoustic_transient_planar_linear_form_2_1 : public VectorFormVol<Scalar>
{
public:
  volume_vector_acoustic_transient_planar_linear_form_2_1(unsigned int i, double ac_rho, double ac_vel);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;

  VectorFormVol<Scalar>* clone() const;

  double ac_rho;
  double ac_vel;

};

template<typename Scalar>
class surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance : public MatrixFormSurf<Scalar>
{
public:
  surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance(unsigned int i, unsigned int j);

  virtual Scalar value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const;
  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const;
  MatrixFormSurf<Scalar>* clone() const;

  double ac_Z0;
};

#pragma endregion

class MyWeakForm : public WeakForm<double>
{
public:
  MyWeakForm(
    Hermes::vector<std::string> acoustic_impedance_markers,
    Hermes::vector<double> acoustic_impedance_values,
    Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns
    ) : WeakForm<double>(2), acoustic_impedance_markers(acoustic_impedance_markers), acoustic_impedance_values(acoustic_impedance_values), prev_slns(prev_slns)
  {
    this->set_ext(prev_slns);

    double acoustic_density = 1.25;
    double acoustic_speed = 343.;
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_1_1<double>(0, 0, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_1_2<double>(0, 1, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_2_1<double>(1, 0, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_2_2<double>(1, 1, acoustic_density, acoustic_speed));

    this->add_vector_form(new volume_vector_acoustic_transient_planar_linear_form_1_2<double>(0, acoustic_density, acoustic_speed));
    this->add_vector_form(new volume_vector_acoustic_transient_planar_linear_form_2_1<double>(1, acoustic_density, acoustic_speed));

    for(int i = 0; i < acoustic_impedance_markers.size(); i++)
    {
      surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<double>* matrix_form = new surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<double>(0, 1);
      matrix_form->set_area(acoustic_impedance_markers[i]);
      matrix_form->ac_Z0 = acoustic_impedance_values[i];
      this->add_matrix_form_surf(matrix_form);
    }
  }

  MyWeakForm(
    std::string acoustic_impedance_marker,
    double acoustic_impedance_value,
    Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns
    ) : WeakForm<double>(2), prev_slns(prev_slns)
  {
    this->set_ext(prev_slns);

    double acoustic_density = 1.25;
    double acoustic_speed = 343.;
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_1_1<double>(0, 0, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_1_2<double>(0, 1, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_2_1<double>(1, 0, acoustic_density, acoustic_speed));
    this->add_matrix_form(new volume_matrix_acoustic_transient_planar_linear_form_2_2<double>(1, 1, acoustic_density, acoustic_speed));

    this->add_vector_form(new volume_vector_acoustic_transient_planar_linear_form_1_2<double>(0, acoustic_density, acoustic_speed));
    this->add_vector_form(new volume_vector_acoustic_transient_planar_linear_form_2_1<double>(1, acoustic_density, acoustic_speed));

    surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<double>* matrix_form = new surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<double>(0, 1);
    matrix_form->set_area(acoustic_impedance_marker);
    matrix_form->ac_Z0 = acoustic_impedance_value;
    this->add_matrix_form_surf(matrix_form);
  }

  Hermes::vector<std::string> acoustic_impedance_markers;
  Hermes::vector<double> acoustic_impedance_values;
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns;
};

class CustomBCValue : public EssentialBoundaryCondition<double>
{
public:
  CustomBCValue(Hermes::vector<std::string> markers, double amplitude = 1., double frequency = 1000.)
    : EssentialBoundaryCondition<double>(markers), amplitude(amplitude), frequency(frequency)
  {
  }

  CustomBCValue(std::string marker, double amplitude = 1., double frequency = 1000.)
    : EssentialBoundaryCondition<double>(marker), amplitude(amplitude), frequency(frequency)
  {
  }

  inline EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return this->amplitude * std::sin(2 * M_PI * this->frequency * this->current_time);
  }

  double amplitude;
  double frequency;
};

class CustomBCDerivative : public EssentialBoundaryCondition<double>
{
public:
  CustomBCDerivative(Hermes::vector<std::string> markers, double amplitude = 1., double frequency = 1000.)
    : EssentialBoundaryCondition<double>(markers), amplitude(amplitude), frequency(frequency)
  {
  }

  CustomBCDerivative(std::string marker, double amplitude = 1., double frequency = 1000.)
    : EssentialBoundaryCondition<double>(marker), amplitude(amplitude), frequency(frequency)
  {
  }

  inline EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return 2 * M_PI * this->frequency * this->amplitude * std::cos(2 * M_PI * this->frequency * this->current_time);
  }

  double amplitude;
  double frequency;
};


template<typename Scalar>
class exact_acoustic_transient_planar_linear_form_1_0_acoustic_pressure : public ExactSolutionScalar<Scalar>
{
public:
  exact_acoustic_transient_planar_linear_form_1_0_acoustic_pressure(MeshSharedPtr mesh, double amplitude, double frequency) : ExactSolutionScalar<Scalar>(mesh), amplitude(amplitude), frequency(frequency)
  {

  }

    Scalar value(double x, double y) const
    {
      return this->amplitude * std::sin(2 * M_PI * this->frequency * this->time);
    }

    void derivatives (double x, double y, Scalar& dx, Scalar& dy) const {};

    Hermes::Ord ord (Hermes::Ord x, Hermes::Ord y) const
    {
        return Hermes::Ord(Hermes::Ord::get_max_order());
    }

    double amplitude;
    double frequency;
    double time;
};


template<typename Scalar>
class exact_acoustic_transient_planar_linear_form_2_0_acoustic_pressure : public ExactSolutionScalar<Scalar>
{
public:
    exact_acoustic_transient_planar_linear_form_2_0_acoustic_pressure(MeshSharedPtr mesh, double amplitude, double frequency) : ExactSolutionScalar<Scalar>(mesh), amplitude(amplitude), frequency(frequency)
    {

    }

    Scalar value(double x, double y) const
    {
       return 2 * M_PI * this->frequency * this->amplitude * std::cos(2 * M_PI * this->frequency * this->time);
    }

    void derivatives (double x, double y, Scalar& dx, Scalar& dy) const {};

    Hermes::Ord ord (Hermes::Ord x, Hermes::Ord y) const
    {
        return Hermes::Ord(Hermes::Ord::get_max_order());
    }

    double amplitude;
    double frequency;
    double time;
};