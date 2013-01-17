#include "parameter/parameter_elemwise_mesh_func.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ParameterElemwiseMeshFunc<Scalar>::ParameterElemwiseMeshFunc(MeshFunction<Scalar>* mesh_function) : ParameterElemwise<Scalar>(), mesh_function(mesh_function)
    {
    }

    template<typename Scalar>
    ParameterElemwiseMeshFunc<Scalar>::ParameterElemwiseMeshFunc() : ParameterElemwise<Scalar>(), mesh_function(NULL)
    {
    }

    template<typename Scalar>
    ParameterElemwiseMeshFunc<Scalar>::~ParameterElemwiseMeshFunc()
    {
    }

    template<typename Scalar>
    ParameterElemwiseValueType ParameterElemwiseMeshFunc<Scalar>::get_type()
    {
      return ParameterElemwiseMeshFunctionValue;
    }

    template<typename Scalar>
    Scalar ParameterElemwiseMeshFunc<Scalar>::get_value(double x, double y)
    {
      return this->mesh_function->get_pt_value()->val[0];
    }

    template<typename Scalar>
    bool ParameterElemwiseMeshFunc<Scalar>::isOkay() const
    {
      if(this->mesh_function == NULL)
        return false;
      return true;
    }
  }
}
