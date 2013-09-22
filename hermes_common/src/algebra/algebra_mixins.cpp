// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file algebra_mixins.cpp
\brief Algebra mixins
*/
#include "algebra_mixins.h"
#include "vector.h"
#include "matrix.h"

namespace Hermes
{
  namespace Algebra
  {
    namespace Mixins
    {
      template<typename Scalar>
      MatrixRhsOutput<Scalar>::MatrixRhsOutput() : output_matrixOn(false), output_matrixIterations(-1), matrixFilename("Matrix_"),
        matrixVarname("A"), matrixFormat(Hermes::Algebra::EXPORT_FORMAT_PLAIN_ASCII), matrix_number_format("%lf"), output_rhsOn(false), output_rhsIterations(-1),
        RhsFilename("Rhs_"), RhsVarname("b"), RhsFormat(Hermes::Algebra::EXPORT_FORMAT_PLAIN_ASCII), rhs_number_format("%lf")
      {
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix, int iteration)
      {
        if (matrix == NULL)
          return;

        char* fileName = new char[this->matrixFilename.length() + 5];

        if(this->output_matrixOn)
        {
          if(this->only_lastMatrixIteration)
            sprintf(fileName, "%s", this->matrixFilename.c_str());
          else if(this->output_matrixIterations == -1 || this->output_matrixIterations >= iteration)
            sprintf(fileName, "%s%i", this->matrixFilename.c_str(), iteration);
          else
          {
            delete [] fileName;
            return;
          }

          matrix->export_to_file(fileName, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
          delete [] fileName;
        }
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::process_matrix_output(Hermes::Algebra::SparseMatrix<Scalar>* matrix)
      {
        if (matrix == NULL)
          return;

        if(this->output_matrixOn)
        {
          char* fileName = new char[this->matrixFilename.length() + 5];
          sprintf(fileName, "%s", this->matrixFilename.c_str());
          matrix->export_to_file(fileName, this->matrixVarname.c_str(), this->matrixFormat, this->matrix_number_format);
          delete [] fileName;
        }
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs, int iteration)
      {
        if (rhs == NULL)
          return;

        char* fileName = new char[this->RhsFilename.length() + 5];

        if(this->output_rhsOn)
        {
          if(this->only_lastRhsIteration)
            sprintf(fileName, "%s", this->RhsFilename.c_str());
          else if(this->output_rhsIterations == -1 || this->output_rhsIterations >= iteration)
            sprintf(fileName, "%s%i", this->RhsFilename.c_str(), iteration);
          else
          {
            delete [] fileName;
            return;
          }

          rhs->export_to_file(fileName, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          delete [] fileName;
        }
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::process_vector_output(Hermes::Algebra::Vector<Scalar>* rhs)
      {
        if (rhs == NULL)
          return;

        if(this->output_rhsOn)
        {
          char* fileName = new char[this->RhsFilename.length() + 5];
          sprintf(fileName, "%s", this->RhsFilename.c_str());
          rhs->export_to_file(fileName, this->RhsVarname.c_str(), this->RhsFormat, this->rhs_number_format);
          delete [] fileName;
        }
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::output_matrix(bool only_last_iteration, int firstIterations)
      {
        output_matrixOn = true;
        this->only_lastMatrixIteration = only_last_iteration;
        this->output_matrixIterations = firstIterations;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_filename(std::string name)
      {
        this->matrixFilename = name;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_print_zero_matrix_entries(bool to_set)
      {
        this->print_matrix_zero_values = to_set;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_varname(std::string name)
      {
        this->matrixVarname = name;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_export_format(Hermes::Algebra::MatrixExportFormat format)
      {
        this->matrixFormat = format;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_matrix_number_format(char* number_format)
      {
        this->matrix_number_format = number_format;
      }

      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::output_rhs(bool only_last_iteration, int firstIterations)
      {
        this->output_rhsOn = true;
        this->only_lastRhsIteration = only_last_iteration;
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
      void MatrixRhsOutput<Scalar>::set_rhs_export_format(Hermes::Algebra::MatrixExportFormat format)
      {
        this->RhsFormat = format;
      }
      template<typename Scalar>
      void MatrixRhsOutput<Scalar>::set_rhs_number_format(char* number_format)
      {
        this->rhs_number_format = number_format;
      }

      template<typename Scalar>
      void MatrixRhsImportExport<Scalar>::export_to_file(std::string filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
      {
        this->export_to_file(filename.c_str(), var_name, fmt, number_format);
      }

      template<typename Scalar>
      void MatrixRhsImportExport<Scalar>::export_to_file(std::string filename, std::string var_name, MatrixExportFormat fmt, char* number_format)
      {
        this->export_to_file(filename.c_str(), var_name.c_str(), fmt, number_format);
      }

      template<typename Scalar>
      void MatrixRhsImportExport<Scalar>::import_from_file(std::string filename, const char *var_name, MatrixExportFormat fmt)
      {
        this->import_from_file(filename.c_str(), var_name, fmt);
      }

      template<typename Scalar>
      void MatrixRhsImportExport<Scalar>::import_from_file(std::string filename, std::string var_name, MatrixExportFormat fmt)
      {
        this->import_from_file(filename.c_str(), var_name.c_str(), fmt);
      }

      template HERMES_API class MatrixRhsOutput<double>;
      template HERMES_API class MatrixRhsOutput<std::complex<double> >;
      template HERMES_API class MatrixRhsImportExport<double>;
      template HERMES_API class MatrixRhsImportExport<std::complex<double> >;
    }
  }
}
