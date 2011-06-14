#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al, HermesFunction<double>* lambda_al,
                        std::string mat_cu, HermesFunction<double>* lambda_cu,
                        HermesFunction<double>* src_term);
};
