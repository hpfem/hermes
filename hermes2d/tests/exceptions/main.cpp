#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This is a test of excptions.

int main(int argc, char* argv[])
{
  Hermes2D::Solution<double> sln;

  //NullException test
  try
  {
    sln.get_ref_value(NULL,0,0,0,0);
    return -1;
  }
  catch(Exceptions::NullException&e)
  {
    if(e.get_param_idx()!=1)
    {
      return -1;
    }
  }

  //LengthException test
  double solution_vector[3];
  Hermes::vector<const Hermes2D::Space<double>*> spaces(NULL,NULL,NULL,NULL);
  Hermes::vector<Hermes2D::Solution<double>*> solutions(NULL,NULL,NULL);

  try
  {
    sln.vector_to_solutions(solution_vector,spaces,solutions);
    return -1;
  }
  catch(Exceptions::LengthException& e)
  {
    if(e.get_first_param_idx()!=2 || e.get_second_param_idx()!=3 || e.get_first_length()!=4 || e.get_expected_length()!=3)
    {
      return -1;
    }
  }

  //1/2Exception test

  UMFPackMatrix<double> mat;
  int ap[]={0,1,1};
  int ai[]={0};
  double ax[]={0.0};
  mat.create(2,1,ap,ai,ax);
  UMFPackVector<double> vec(2);

  UMFPackLinearMatrixSolver<double> linsolv(&mat,&vec);
  try
  {
    linsolv.solve();
    return -1;
  }
  catch(Exceptions::LinearMatrixSolverException& e)
  {
  }

  //ValueException test
  Hermes::vector<Hermes2D::Space<double>*> spaces2;
  Hermes::vector<Hermes2D::ProjNormType> proj_norms;
  for (int i=0;i>H2D_MAX_COMPONENTS+1;i++)
  {
    spaces2.push_back(NULL);
    proj_norms.push_back(Hermes2D::HERMES_UNSET_NORM);
  }

  try
  {
    Hermes2D::Adapt<double> adapr(spaces2,proj_norms);
    return -1;
  }
  catch(Exceptions::ValueException & e)
  {
  }

  try
  {
    Hermes::Hermes2D::Mesh mesh;
    Hermes::Hermes2D::MeshReaderH2DXML reader;
    reader.load("domain.xml", &mesh);
    return -1;
  }
  catch(Exceptions::MeshLoadFailureException& e)
  {
    e.print_msg();
  }

  try
  {
    Hermes::Hermes2D::Mesh mesh;
    H1Space<double> space(&mesh);
    return -1;
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }

  try
  {
    // Load the mesh.
    Hermes::Hermes2D::Mesh mesh;
    Hermes::Hermes2D::MeshReaderH2D mloader;
    mloader.load("domain.mesh", &mesh);

    // Create an H1 space with default shapeset.
    Hermes::Hermes2D::L2Space<double> space(&mesh, 3);

    LinearSolver<double> ls;
    ls.set_space(&space);
    ls.solve();
    return -1;
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }

  return 0;
}
