// This file is part of Hermes2D.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, see <http://www.gnu.prg/licenses/>.
#include "mixins2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Mixins
    {
      template<typename Scalar>
      const Space<Scalar>* SettableSpaces<Scalar>::get_space(int n) const
      {
        return this->get_spaces()[n];
      }

      void StateQueryable::check() const
      {
        if(!this->isOkay())
          {
            std::stringstream ss;
            ss << "The instance of " << this->getClassName() << " is not OK.";
            throw Hermes::Exceptions::Exception(ss.str().c_str());
          }
      }

      XMLParsing::XMLParsing() : validate(true)
      {
      }

      void XMLParsing::set_validation(bool to_set)
      {
        this->validate = to_set;
      }

      template<typename Scalar>
      MatrixRhsOutput<Scalar>::MatrixRhsOutput() : output_matrixOn(false), output_matrixIterations(-1), matrixFilename("Matrix_"),
        matrixVarname("A"), matrixFormat(Hermes::Algebra::DF_MATLAB_SPARSE), matrix_number_format("%lf"), output_rhsOn(false), output_rhsIterations(-1),
        RhsFilename("Rhs_"), RhsVarname("b"), RhsFormat(Hermes::Algebra::DF_MATLAB_SPARSE), rhs_number_format("%lf")
      {
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::output_matrix(int firstIterations)
      {
        output_matrixOn = true;
        this->output_matrixIterations = firstIterations;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_filename(std::string name)
      {
        this->matrixFilename = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_varname(std::string name)
      {
        this->matrixVarname = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_E_matrix_dump_format(EMatrixDumpFormat format)
      {
        this->matrixFormat = format;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_number_format(char* number_format)
      {
        this->matrix_number_format = number_format;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::output_rhs(int firstIterations)
      {
        this->output_rhsOn = true;
        this->output_rhsIterations = firstIterations;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_rhs_filename(std::string name)
      {
        this->RhsFilename = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_rhs_varname(std::string name)
      {
        this->RhsVarname = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_rhs_E_matrix_dump_format(EMatrixDumpFormat format)
      {
        this->RhsFormat = format;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_rhs_number_format(char* number_format)
      {
        this->rhs_number_format = number_format;
      }

      template HERMES_API class SettableSpaces<double>;
      template HERMES_API class SettableSpaces<std::complex<double> >;
      template HERMES_API class MatrixRhsOutput<double>;
      template HERMES_API class MatrixRhsOutput<std::complex<double> >;
    }
  }
}
