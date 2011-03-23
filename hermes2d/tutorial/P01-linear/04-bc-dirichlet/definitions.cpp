#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double const_f) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormVolConst(0, 0));
    add_vector_form(new DefaultVectorFormVolConst(0, const_f));
  };
};

class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker, double const_f) 
           : EssentialBC(marker), const_f(const_f) { }

  ~EssentialBCNonConst() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return (-const_f/4.0)*(x*x + y*y);
  }

  // Member.
  double const_f;
};
