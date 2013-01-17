#include "elemwise_parameter/elemwise_parameter_mesh_func.h"
#include "../forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ElemwiseParameterMeshFunc<Scalar>::ElemwiseParameterMeshFunc(MeshFunction<Scalar>* mesh_function) : ElemwiseParameterFunc<Scalar>(), mesh_function(mesh_function)
    {
    }

    template<typename Scalar>
    ElemwiseParameterMeshFunc<Scalar>::ElemwiseParameterMeshFunc() : ElemwiseParameterFunc<Scalar>(), mesh_function(NULL)
    {
    }

    template<typename Scalar>
    ElemwiseParameterMeshFunc<Scalar>::~ElemwiseParameterMeshFunc()
    {
    }

    template<typename Scalar>
    Scalar ElemwiseParameterMeshFunc<Scalar>::get_value(Element* e)
    {
      double x, y;
      e->get_center(x, y);
      return this->mesh_function->get_pt_value(x, y)->val[0];
    }

    template<typename Scalar>
    bool ElemwiseParameterMeshFunc<Scalar>::isOkay() const
    {
      if(this->mesh_function == NULL)
        return false;
      return true;
    }

    template<typename Scalar>
    std::string ElemwiseParameterMeshFunc<Scalar>::getClassName() const
    {
      return "ElemwiseParameterMeshFunc";
    }

    template class HERMES_API ElemwiseParameterMeshFunc<double>;
    template class HERMES_API ElemwiseParameterMeshFunc<std::complex<double> >;
  }
}
