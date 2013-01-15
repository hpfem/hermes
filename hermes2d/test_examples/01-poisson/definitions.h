#include "hermes2d.h"

/* Weak forms */
using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoisson();
};

class CustomWeakFormPoissonCacheCalculation : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoissonCacheCalculation();
};
