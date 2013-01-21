#include "elemwise_parameter/elemwise_parameter_nonlinear.h"
#include "forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ElemwiseParameterNonlinear<Scalar>::ElemwiseParameterNonlinear(ElemwiseParameterNonlinearValueType value_type) : ElemwiseParameter<Scalar>(), value_type(value_type)
    {
    }

    template<typename Scalar>
    ElemwiseParameterNonlinear<Scalar>::~ElemwiseParameterNonlinear()
    {
    }

    template<typename Scalar>
    ElemwiseParameterType ElemwiseParameterNonlinear<Scalar>::get_type()
    {
      return ElemwiseParameterTypeNonlinear;
    }

    template<typename Scalar>
    ElemwiseParameterNonlinearHermesFunc<Scalar>::ElemwiseParameterNonlinearHermesFunc(Hermes1DFunction<Scalar>* function, ElemwiseParameterNonlinearValueType value_type) : ElemwiseParameterNonlinear<Scalar>(value_type), function(function)
    {
    }

    template<typename Scalar>
    ElemwiseParameterNonlinearHermesFunc<Scalar>::ElemwiseParameterNonlinearHermesFunc(ElemwiseParameterNonlinearValueType value_type) : ElemwiseParameterNonlinear<Scalar>(value_type)
    {
    }

    template<typename Scalar>
    ElemwiseParameterNonlinearHermesFunc<Scalar>::~ElemwiseParameterNonlinearHermesFunc()
    {
    }

    template<typename Scalar>
    Scalar ElemwiseParameterNonlinearHermesFunc<Scalar>::get_value(int np, Func<Scalar>* u_ext, Geom<double>* geometry)
    {
      Scalar result = 0.0;
      for(int i = 0; i < np; i++)
        result += u_ext->val[i];
      result /= np;
      
      // Pass the value.
      if(this->value_type == ElemwiseParameterNonlinearValue)
        return this->function->value(result);
      else
        return this->function->derivative(result);
    }

    template<typename Scalar>
    bool ElemwiseParameterNonlinearHermesFunc<Scalar>::isOkay() const
    {
      if(this->function == NULL)
        return false;
      return true;
    }

    template<typename Scalar>
    std::string ElemwiseParameterNonlinearHermesFunc<Scalar>::getClassName() const
    {
      return "ElemwiseParameterNonlinearHermesFunc";
    }

    ElemwiseParameterSpline::ElemwiseParameterSpline(ElemwiseParameterNonlinearValueType value_type) : ElemwiseParameterNonlinearHermesFunc<double>(value_type)
    {
    }

    ElemwiseParameterSpline::ElemwiseParameterSpline(CubicSpline* spline, ElemwiseParameterNonlinearValueType value_type) : ElemwiseParameterNonlinearHermesFunc<double>(spline, value_type)
    {
    }

    ElemwiseParameterSpline::~ElemwiseParameterSpline()
    {
    }

    template class HERMES_API ElemwiseParameterNonlinear<double>;
    template class HERMES_API ElemwiseParameterNonlinear<std::complex<double> >;
    template class HERMES_API ElemwiseParameterNonlinearHermesFunc<double>;
    template class HERMES_API ElemwiseParameterNonlinearHermesFunc<std::complex<double> >;
  }
}
