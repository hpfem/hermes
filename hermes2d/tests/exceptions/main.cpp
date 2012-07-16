#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace Hermes;

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
    if(e.getParamIdx()!=1)
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
    if(e.getFirstParamIdx()!=2 || e.getSecondParamIdx()!=3 || e.getFirstLength()!=4 || e.getExpectedLength()!=3)
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
  }
  catch(Exceptions::MeshLoadFailureException& e)
  {
    e.printMsg();
  }

  return 0;
}