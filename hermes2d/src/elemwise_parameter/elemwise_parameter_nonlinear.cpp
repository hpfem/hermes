#include "elemwise_parameter/elemwise_parameter_nonlinear.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ElemwiseParameterNonlinear<Scalar>::ElemwiseParameterNonlinear() : ElemwiseParameter<Scalar>()
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

    ElemwiseParameterSpline::ElemwiseParameterSpline() : ElemwiseParameterNonlinear(), spline(NULL)
    {
    }

    ElemwiseParameterSpline::ElemwiseParameterSpline(CubicSpline* spline) : ElemwiseParameterNonlinear(), spline(spline)
    {
    }

    ElemwiseParameterSpline::~ElemwiseParameterSpline()
    {
    }

    double ElemwiseParameterSpline::get_value(double u_ext_value)
    {
      int interval;
      if(!spline->find_interval(u_ext_value, interval))
        throw Hermes::Exceptions::Exception("Asked for a spline value out of the interval");

      return spline->get_value_from_interval(u_ext_value, interval);
    }

    bool ElemwiseParameterSpline::isOkay() const
    {
      if(this->spline == NULL)
        return false;
      return true;
    }

    std::string ElemwiseParameterSpline::getClassName() const
    {
      return "ElemwiseParameterSpline";
    }

    template class HERMES_API ElemwiseParameterNonlinear<double>;
    template class HERMES_API ElemwiseParameterNonlinear<std::complex<double> >;
  }
}
