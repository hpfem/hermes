#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

template<typename Scalar>
class ExactSolutionScalarXY : public ExactSolutionScalar<Scalar>
{
public:
  ExactSolutionScalarXY(MeshSharedPtr mesh) : ExactSolutionScalar<Scalar>(mesh) {}
  virtual ~ExactSolutionScalarXY() {};

  /// Function returning the value.
  virtual Scalar value (double x, double y) const
  {
    return x*x - y*y;
  }

  virtual Ord ord(double x, double y) const
  {
    return Hermes::Ord(5);
  }

  MeshFunction<Scalar>* clone() const
  {
    return new ExactSolutionScalarXY<Scalar>(this->mesh);
  }

  /// Function returning the derivatives.
  virtual void derivatives (double x, double y, Scalar& dx, Scalar& dy) const
  {
    dx = 2.0*x;
    dy = -2.0*y;
  }
};

class MagneticVolumetricIntegralEggShellCalculator : public Hermes::Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>
{
public:
  MagneticVolumetricIntegralEggShellCalculator(Hermes::Hermes2D::MeshFunctionSharedPtr<double> source_function, int number_of_integrals)
    : Hermes::Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MagneticVolumetricIntegralEggShellCalculator(Hermes::vector<Hermes::Hermes2D::MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals)
    : Hermes::Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Hermes::Hermes2D::Func<double> **fns, Hermes::Hermes2D::Geom<double> *e, double* result)
  {
    double *x = e->x;
    double *y = e->y;

    // functions
    double **value = new double*[source_functions.size()];
    double **dudx = new double*[source_functions.size()];
    double **dudy = new double*[source_functions.size()];

    for (int i = 0; i < source_functions.size(); i++)
    {
      value[i] = fns[i]->val;
      dudx[i] = fns[i]->dx;
      dudy[i] = fns[i]->dy;
    }

    // expressions
    for (int i = 0; i < n; i++)
    {
      result[0] += wt[i];
      result[1] += wt[i] * (y[i]*((dudx[source_functions.size() - 1][i]*(-1e6*(dudx[0][i]*dudx[0][i]+dudy[0][i]*dudy[0][i])+1e6*(dudx[0][i]*dudx[0][i])))+(dudy[source_functions.size() - 1][i]*1e6*dudx[0][i]*dudy[0][i]))-x[i]*((dudy[source_functions.size() - 1][i]*(-1e6*(dudx[0][i]*dudx[0][i]+dudy[0][i]*dudy[0][i])+1e6*(dudy[0][i]*dudy[0][i])))+(dudx[source_functions.size() - 1][i]*1e6*dudx[0][i]*dudy[0][i])));
    }

    delete [] value;
    delete [] dudx;
    delete [] dudy;
  }

  virtual void order(Hermes::Hermes2D::Func<Hermes::Ord> **fns, Hermes::Ord* result)
  {
    result[0] = Hermes::Ord(20);
    result[1] = Hermes::Ord(20);
  }

};

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  // Hermes::Hermes2D::MeshReaderH2DBSON mloader_bson;
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);
  mloader.load("acoustic.msh", meshes);

  // mloader_bson.save("bson_mesh", meshes);
  // mesh->free();
  // mloader_bson.load("bson_mesh", meshes);


  mesh.get()->refine_towards_boundary("2", 1);
  mesh.get()->refine_towards_boundary("1", 1);
  mesh.get()->refine_towards_boundary("25", 2);

  mesh.get()->refine_all_elements(1);

  for(int asdf = 0; asdf < 20; asdf++)
  {
    Hermes::Hermes2D::MeshSharedPtr eggShellMesh = Hermes::Hermes2D::EggShell::get_egg_shell(mesh, "3", 3);

    MeshView m;
    m.show(eggShellMesh);

    Hermes::Hermes2D::MeshFunctionSharedPtr<double> b;
    Hermes::Hermes2D::MeshFunctionSharedPtr<double> a(new Hermes::Hermes2D::ExactSolutionEggShell(eggShellMesh, 3));

    ScalarView s;
    s.show(a);

    b = a;

    Hermes::Hermes2D::MeshFunctionSharedPtr<double> sln(new ExactSolutionScalarXY<double>(mesh));

    Hermes::vector<Hermes::Hermes2D::MeshFunctionSharedPtr<double> > slns;
    slns.push_back(sln);
    slns.push_back(a);

    MagneticVolumetricIntegralEggShellCalculator calcEggShell(slns, 2);
    double *valuesEggShell = calcEggShell.calculate(Hermes::vector<std::string>("0", "1", "2", "4", "5"));
    std::cout << valuesEggShell[0] << std::endl << valuesEggShell[1] << std::endl;

    MeshView mview("Hello world!", new WinGeom(0, 0, 350, 350));
    mview.show(mesh);
  }
  return 0;
}
