#include "weakform/weakform.h"
#include "weakform_library/maxwell.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "function/function.h"

using namespace WeakFormsMaxwell::VolumetricMatrixForms;
using namespace WeakFormsMaxwell::VolumetricVectorForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm
{
public:
  CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2, 
                               CubicSpline* mu_inv_iron, std::string material_air,
                               std::string material_copper, double mu_vacuum, 
                               double current_density, int order_inc = 3) 
    : WeakForm(1) {

    /* AXISYM LINEAR
    // Jacobian.
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_air, 1.0, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_copper, 1.0, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_iron_1, 1/300.,  
                                                               HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_iron_2, 1/300., 
                                                               HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));

    // Residual.
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_air, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_copper, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_iron_1, 1/300., HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_iron_2, 1/300., HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultVectorFormConst(0, material_copper, -current_density * mu_vacuum, HERMES_PLANAR));
    */

    // AXISYM NONLINEAR
    // Jacobian.
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_air, 1.0, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_copper, 1.0, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_1, mu_inv_iron, 1.0, 
                                                               HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_2, mu_inv_iron, 1.0, 
                                                               HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));

    // Residual.
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_air, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_copper, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_1, mu_inv_iron, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_2, mu_inv_iron, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultVectorFormConst(0, material_copper, -current_density * mu_vacuum, HERMES_PLANAR));
   };
};

class HERMES_API FilterVectorPotencial : public MagFilter
{
public:
  FilterVectorPotencial(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items) : MagFilter(solutions, items) {};

protected:
  void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result, Geom<double> *e)
  {
    for (int i = 0; i < n; i++){
      result[i] = 0;
      for(unsigned int j = 0; j < values.size(); j++)
        result[i] += sqr(values[j][i]);

      result[i] = sqrt(result[i]) * e->x[i];
    }
  }
};

class HERMES_API FilterFluxDensity : public Filter
{
public:
  FilterFluxDensity(Hermes::vector<MeshFunction*> solutions)
           : Filter(solutions) {};

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL) {
    error("Not implemented yet"); return 0;
  }

protected:
  void precalculate(int order, int mask)
  {
    Quad2D* quad = quads[cur_quad];
    int np = quad->get_num_points(order);
    Node* node = new_node(H2D_FN_DEFAULT, np);

    sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
    sln[1]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);

    scalar *dudx1, *dudy1, *dudx2, *dudy2;
    sln[0]->get_dx_dy_values(dudx1, dudy1);
    sln[1]->get_dx_dy_values(dudx2, dudy2);
    scalar *uval1 = sln[0]->get_fn_values();
    scalar *uval2 = sln[1]->get_fn_values();
    update_refmap();
    double *x = refmap->get_phys_x(order);

    for (int i = 0; i < np; i++)
    {
      node->values[0][0][i] = sqrt(sqr(dudy1[i]) + sqr(dudy2[i]) +
                                   sqr(dudx1[i] + ((x[i] > 1e-10) ? uval1[i] / x[i] : 0.0)) +
                                   sqr(dudx2[i] + ((x[i] > 1e-10) ? uval2[i] / x[i] : 0.0)));
    }

    if(nodes->present(order)) {
      assert(nodes->get(order) == cur_node);
      ::free(nodes->get(order));
    }
    nodes->add(node, order);
    cur_node = node;
  }
};
