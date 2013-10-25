#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  Hermes::Hermes2D::MeshReaderH2DBSON mloader_bson;
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh);
  mloader.load("acoustic.msh", meshes);
  
  mloader_bson.save("bson_mesh", meshes);
  mesh->free();
  mloader_bson.load("bson_mesh", meshes);

  Hermes::Hermes2D::MeshSharedPtr eggShellMesh = Hermes::Hermes2D::EggShell::get_egg_shell(mesh, "3", 2);
  Hermes::Hermes2D::MeshFunctionSharedPtr<double> b;
  Hermes::Hermes2D::MeshFunctionSharedPtr<double> a(new Hermes::Hermes2D::ExactSolutionEggShell(eggShellMesh, 3));
  b = a;

  return 0;
}