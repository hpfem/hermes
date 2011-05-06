#include "hermes2d.h"

using namespace WeakFormsH1;

/* Weak forms */

class WeakFormEigenLeft : public WeakForm
{
public:
  WeakFormEigenLeft() : WeakForm(1) {
    add_matrix_form(new DefaultJacobianDiffusion(0, 0));
  };
};

class WeakFormEigenRight : public WeakForm
{
public:
  WeakFormEigenRight() : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormVol(0, 0));
  };
};
