#include "elemwise_parameter/elemwise_parameter_nonlinear.h"

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
    Scalar ElemwiseParameterNonlinearHermesFunc<Scalar>::get_value(Solution<Scalar>* u_ext, Element* e)
    {
      if(u_ext->num_components > 1)
        throw Hermes::Exceptions::Exception("ElemwiseParameterNonlinearHermesFunc does not yet work with vector shapesets.");

      // Get center of the element.
      double x, y;
      e->get_center(x, y);

      // Get the value at x, y.
      int o = u_ext->elem_orders[e->id];
      Scalar* mono = u_ext->dxdy_coeffs[0][0];
      Scalar result = 0.0;
      int k = 0;
      for (int i = 0; i <= o; i++)
      {
        Scalar row = mono[k++];
        for (int j = 0; j < (e->is_quad() ? o : i); j++)
          row = row * x + mono[k++];
        result = result * y + row;
      }

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
