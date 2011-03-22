#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

class WeakFormNewton : public WeakFormLaplace
{
public:
  // Problem parameters.
  WeakFormNewton(double h, double T0, std::string natural_bc_bnd_part) : WeakFormLaplace(1)
  {
    add_matrix_form(new WeakFormLaplace::MatrixFormVol(0, 0));
    add_matrix_form_surf(new WeakFormLaplace::MatrixFormSurfNewton(0, 0, natural_bc_bnd_part, h));
    add_vector_form_surf(new WeakFormLaplace::VectorFormSurfNewton(0, natural_bc_bnd_part, h * T0));
  };
};

