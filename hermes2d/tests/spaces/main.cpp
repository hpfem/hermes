#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace Hermes::Hermes2D;

// This is a test of spaces in Hermes2D.

class CustomExactSol : public ExactSolutionScalar<double>
{
public:
  CustomExactSol(Mesh* mesh_lshape) : ExactSolutionScalar<double>(mesh_lshape) {};

  double value (double x, double y) const { return std::pow(x + y, 2.0); };

  void derivatives (double x, double y, double& dx, double& dy) const {};

  Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const { return Hermes::Ord(1); }
};

class BCNonConst : public DefaultEssentialBCNonConst<double>
{
public:
  BCNonConst(std::string marker, ExactSolutionScalar<double>* exact_solution) : DefaultEssentialBCNonConst<double>(marker, exact_solution) {};

  double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return exact_solution->value(x, y);
  }

  EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }
};

int main(int argc, char* argv[])
{
  // Tested entities.
  Hermes::vector<int> dofs;
  Hermes::vector<double*> bc_values;

  // Known values (tested Element IDs etc.).
  Hermes::vector<Hermes::vector<int> > element_ids_boundary(
    Hermes::vector<int>(4, 21, 101),
    Hermes::vector<int>(13, 46, 170),
    Hermes::vector<int>(17, 72, 293),
    Hermes::vector<int>(19, 71, 291));

  // Load the mesh_lshape.
  Mesh mesh_lshape, mesh_square;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh_lshape);
  mloader.load("square.mesh", &mesh_square);

  // Basic DOFs assignment and BC values testing.
  for(unsigned int refinement_i = 0; refinement_i < 3; refinement_i++)
  {
    // Mesh refinement.
    mesh_lshape.refine_all_elements();
    mesh_square.refine_all_elements();

    // H1Spaces from order 1 for various boundary conditions:
    CustomExactSol exact(&mesh_lshape);
    Hermes::vector<BCNonConst> bc_essential(
      BCNonConst("Bottom", &exact),
      BCNonConst("Outer", &exact),
      BCNonConst("Inner", &exact),
      BCNonConst("Left", &exact));

    for(unsigned int poly_order_i = 1; poly_order_i < 10; poly_order_i++)
      for(unsigned int bc_i = 0; bc_i < 4; bc_i++)
      {
        // Create a space and record its DOFs.
        EssentialBCs<double> bcs(&(bc_essential[bc_i]));
        H1Space<double>* space = new H1Space<double>(&mesh_lshape, &bcs, poly_order_i);

        space->save("space.xml");

        space->load("space.xml", &mesh_lshape, &bcs);

        dofs.push_back(space->get_num_dofs());

        SurfPos surf_pos;
        surf_pos.base = mesh_lshape.get_element_fast(element_ids_boundary[bc_i][refinement_i]);
        if(bc_i == 0) {
          surf_pos.surf_num = 0;
          surf_pos.marker = mesh_lshape.get_boundary_markers_conversion().get_internal_marker("Bottom").marker;
        }
        if(bc_i == 1) {
          // This would be for NURBS boundaries, there is a separate test for this.
          continue;
        }
        if(bc_i == 2) {
          surf_pos.surf_num = 0;
          surf_pos.marker = mesh_lshape.get_boundary_markers_conversion().get_internal_marker("Inner").marker;
        }
        if(bc_i == 3) {
          surf_pos.surf_num = 3;
          surf_pos.marker = mesh_lshape.get_boundary_markers_conversion().get_internal_marker("Left").marker;
        }
        surf_pos.lo = surf_pos.hi = 0;

        bc_values.push_back(space->get_bc_projection(&surf_pos, 3));
      }

      // L2Spaces from order 0.
      for(unsigned int poly_order_i = 0; poly_order_i < 10; poly_order_i++)
      {
        // Create a space and record its DOFs.
        L2Space<double> space(&mesh_lshape, poly_order_i);
        dofs.push_back(space.get_num_dofs());
      }

      // HCurl & HDiv spaces from order 1.
      for(unsigned int poly_order_i = 1; poly_order_i < 10; poly_order_i++)
      {
        // Create a space and record its DOFs.
        HcurlSpace<double> space1(&mesh_lshape, poly_order_i);
        dofs.push_back(space1.get_num_dofs());

        // Create a space and record its DOFs.
        /// \todo HDiv shapeset on triangles?
        HdivSpace<double> space2(&mesh_square, poly_order_i);
        dofs.push_back(space2.get_num_dofs());
      }
  }

  /// \todo Adjusting element orders testing.

  /// \todo Create refined spaces testing.

  /// \todo Space copying testing.

  /// \todo AsmList generating testing.

  /// \todo HCurl & HDiv spaces boundary condition testing.

  /* This is how the data were generated
  std::ofstream dofs_out("test_dofs.txt");
  for(unsigned int dofs_out_i = 0; dofs_out_i < dofs.size(); dofs_out_i++)
    dofs_out << dofs[dofs_out_i] << std::endl;
  dofs_out.close();

  std::ofstream bc_values_out("test_bc_values.txt");
  for(unsigned int bc_values_out_i = 0; bc_values_out_i < bc_values.size(); bc_values_out_i++)
    for(unsigned int value_i = 0; value_i < 4; value_i++)
      bc_values_out << bc_values[bc_values_out_i][value_i] << std::endl;
  bc_values_out.close();
  */

  // Test itself.
  std::ifstream dofs_in("test_dofs.txt");
  int dofs_test;
  for(unsigned int dofs_in_i = 0; dofs_in_i < dofs.size(); dofs_in_i++)
  {
    dofs_in >> dofs_test;
    if(dofs_test != dofs[dofs_in_i])
    {
      return -1;
    }
  }
  dofs_in.close();

  std::ifstream bc_values_in("test_bc_values.txt");
  double bc_values_test;
  for(unsigned int bc_values_in_i = 0; bc_values_in_i < bc_values.size(); bc_values_in_i++)
    for(unsigned int value_i = 0; value_i < 4; value_i++)
    {
      bc_values_in >> bc_values_test;
      if(std::abs(bc_values[bc_values_in_i][value_i] - bc_values_test) > 1E-8)
      {
        return -1;
      }
    }
  bc_values_in.close();

  return 0;
}