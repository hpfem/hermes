#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al, Hermes::Hermes2D::HermesFunction<double>* lambda_al,
                        std::string mat_cu, Hermes::Hermes2D::HermesFunction<double>* lambda_cu,
                        Hermes::Hermes2D::HermesFunction<double>* src_term);
};
