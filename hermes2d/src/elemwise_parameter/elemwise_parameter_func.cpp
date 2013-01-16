#include "elemwise_parameter/elemwise_parameter_func.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ElemwiseParameterFunc<Scalar>::ElemwiseParameterFunc() : ElemwiseParameter<Scalar>()
    {
    }

    template<typename Scalar>
    ElemwiseParameterFunc<Scalar>::~ElemwiseParameterFunc()
    {
    }

    template<typename Scalar>
    ElemwiseParameterType ElemwiseParameterFunc<Scalar>::get_type()
    {
      return ElemwiseParameterTypeFunc;
    }

    template class HERMES_API ElemwiseParameterFunc<double>;
    template class HERMES_API ElemwiseParameterFunc<std::complex<double> >;
  }
}
