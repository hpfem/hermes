#include "elemwise_parameter/elemwise_parameter.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ElemwiseParameter<Scalar>::ElemwiseParameter()
    {
    }

    template<typename Scalar>
    ElemwiseParameter<Scalar>::~ElemwiseParameter()
    {
    }

    template class HERMES_API ElemwiseParameter<double>;
    template class HERMES_API ElemwiseParameter<std::complex<double> >;
  }
}
