#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(HERMES_ANY, 10.);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space.
  SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, nullptr, 2));
  space->set_element_order(0, 1);
  space->set_element_order(1, 2);
  space->assign_dofs();
  BaseView<double> b("Coarse space");
  b.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  b.show(space);

  MeshSharedPtr meshf(new Mesh);
  meshf->copy(mesh);
  meshf->refine_element(meshf->get_element(0), 0);
  SpaceSharedPtr<double> spacef(new H1Space<double>(meshf, 1));
  spacef->set_element_order(1, 2);
  spacef->set_element_order(2, 1);
  spacef->set_element_order(3, 1);
  spacef->set_element_order(4, 1);
  spacef->set_element_order(5, 1);
  spacef->assign_dofs();

  BaseView<double> bf("Fine space");
  bf.get_linearizer()->set_criterion(LinearizerCriterionFixed(3));
  bf.show(spacef);

  View::wait();
  return 0;
}