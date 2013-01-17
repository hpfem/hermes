#include "parameter/parameter_elemwise_func.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ParameterElemwiseFunc<Scalar>::ParameterElemwiseFunc() : ParameterElemwise<Scalar>()
    {
    }

    template<typename Scalar>
    ParameterElemwiseFunc<Scalar>::~ParameterElemwiseFunc()
    {
    }

    template<typename Scalar>
    ParameterElemwiseValueType ParameterElemwiseFunc<Scalar>::get_type()
    {
      return ParameterElemwiseFunctionValue;
    }
  }
}
