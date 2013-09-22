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
/*! \file vector.cpp
\brief Basic vector classes and operations.
*/
#include "common.h"
#include "matrix.h"
#include "callstack.h"

#include "solvers/linear_matrix_solver.h"
#include "solvers/umfpack_solver.h"
#include "solvers/superlu_solver.h"
#include "solvers/amesos_solver.h"
#include "solvers/petsc_solver.h"
#include "solvers/mumps_solver.h"
#include "solvers/aztecoo_solver.h"
#include "solvers/paralution_solver.h"
#include "qsort.h"
#include "api.h"

namespace Hermes
{
  namespace Algebra
  {
    template<typename Scalar>
    Vector<Scalar>::Vector() : size(0)
    {
    }

    template<typename Scalar>
    Vector<Scalar>::Vector(unsigned int size) : size(size)
    {
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->set(i, vec->get(i));
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::set_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->set(i, vec[i]);
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::add_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i));
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* Vector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i]);
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format)
    {
      if(!v)
        throw Exceptions::MethodNotOverridenException("Vector<Scalar>::export_to_file");

      switch (fmt)
      {
      case EXPORT_FORMAT_MATRIX_MARKET:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate real\n");
          else
            fprintf(file, "%%%%Matrix<Scalar>Market matrix coordinate complex\n");

          fprintf(file, "%d 1 %d\n", this->size, this->size);

          for (unsigned int j = 0; j < this->size; j++)
          {
            Hermes::Helpers::fprint_coordinate_num(file, j + 1, 1, v[j], number_format);
            fprintf(file, "\n");
          }

          fclose(file);
        }
        break;

      case EXPORT_FORMAT_MATLAB_MATIO:
        {
#ifdef WITH_MATIO
          size_t dims[2];
          dims[0] = this->size;
          dims[1] = 1;

          mat_t *mat = Mat_CreateVer(filename, "", MAT_FT_MAT5);
          matvar_t *matvar;

          // For complex.
          double* v_re = NULL;
          double* v_im = NULL;

          void* data;
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
          {
            data = v;
            matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA);
          }
          else
          {
            v_re = new double[this->size];
            v_im = new double[this->size];
            struct mat_complex_split_t z = {v_re, v_im};

            for(int i = 0; i < this->size; i++)
            {
              v_re[i] = ((std::complex<double>)(this->v[i])).real();
              v_im[i] = ((std::complex<double>)(this->v[i])).imag();
              data = &z;
            }
            matvar = Mat_VarCreate(var_name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
          } 

          if (matvar)
          {
            Mat_VarWrite(mat, matvar, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
          }

          if(v_re)
            delete [] v_re;
          if(v_im)
            delete [] v_im;
          Mat_Close(mat);

          if(!matvar)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
#else
          throw Exceptions::Exception("MATIO not included.");
#endif
        }
        break;

      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          FILE* file = fopen(filename, "w");
          if(!file)
            throw Exceptions::IOException(Exceptions::IOException::Write, filename);
          for (unsigned int i = 0; i < this->size; i++)
          {
            Hermes::Helpers::fprint_num(file, v[i], number_format);
            fprintf(file, "\n");
          }
          fclose(file);
        }
      }
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt)
    {
      switch (fmt)
      {
      case EXPORT_FORMAT_PLAIN_ASCII:
        {
          std::vector<Scalar> data;
          std::ifstream input (filename);
          if(input.bad())
            throw Exceptions::IOException(Exceptions::IOException::Read, filename);
          std::string lineData;

          while(getline(input, lineData))
          {
            Scalar d;
            std::stringstream lineStream(lineData);
            lineStream >> d;
            data.push_back(d);
          }

          this->alloc(data.size());
          memcpy(this->v, &data[0], sizeof(Scalar)*data.size());
        }
        break;
      case EXPORT_FORMAT_MATLAB_MATIO:
#ifdef WITH_MATIO
        mat_t    *matfp;
        matvar_t *matvar;

        matfp = Mat_Open(filename,MAT_ACC_RDONLY);

        if (!matfp )
        {
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
          return;
        }

        matvar = Mat_VarRead(matfp, var_name);
        if (matvar)
        {
          this->alloc(matvar->dims[0]);
          if(Hermes::Helpers::TypeIsReal<Scalar>::value)
            memcpy(this->v, matvar->data, sizeof(Scalar)*this->size);
          else
          {
            std::complex<double>* complex_data = new std::complex<double>[this->size];
            double* real_array = (double*)((mat_complex_split_t*)matvar->data)->Re;
            double* imag_array = (double*)((mat_complex_split_t*)matvar->data)->Im;
            for(int i = 0; i < this->size; i++)
              complex_data[i] = std::complex<double>(real_array[i], imag_array[i]);
            memcpy(this->v, complex_data, sizeof(Scalar)*this->size);
            delete [] complex_data;
          }
        }

        Mat_Close(matfp);
        if(!matvar)
          throw Exceptions::IOException(Exceptions::IOException::Read, filename);
#else
        throw Exceptions::Exception("MATIO not included.");
#endif
        break;
      case EXPORT_FORMAT_MATRIX_MARKET:
        throw Hermes::Exceptions::MethodNotImplementedException("SimpleVector<Scalar>::import_from_file - Matrix Market");
      }

    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector() : Vector<Scalar>(), v(NULL)
    {
    }

    template<typename Scalar>
    SimpleVector<Scalar>::SimpleVector(unsigned int size) : Vector<Scalar>(size), v(NULL)
    {
      if(this->size == 0)
        throw Exceptions::ValueException("size", this->size, 1);
      this->alloc(this->size);
    }

    template<typename Scalar>
    SimpleVector<Scalar>::~SimpleVector()
    {
      free();
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::set_vector(Hermes::Algebra::Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      SimpleVector<Scalar>* simple_vec = (SimpleVector<Scalar>*)vec;
      if(simple_vec)
        memcpy(this->v, simple_vec->v, sizeof(Scalar)*this->size);
      else
        Vector<Scalar>::set_vector(vec);
      return this;
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::set_vector(Scalar* vec)
    {
      memcpy(this->v, vec, sizeof(Scalar)*this->size);
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      this->v = new Scalar[n];
      zero();
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++)
        v[i] *= -1.;
      return this;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::zero()
    {
      memset(this->v, 0, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::free()
    {
      if (this->v)
        delete [] this->v;
      this->v = NULL;
      this->size = 0;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      this->v[idx] = y;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add(unsigned int idx, Scalar y)
    {
#pragma omp critical (SimpleVector_add)
      this->v[idx] += y;
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        this->v[idx[i]] += y[i];
    }

    template<typename Scalar>
    Scalar SimpleVector<Scalar>::get(unsigned int idx) const
    {
      return this->v[idx];
    }

    template<typename Scalar>
    void SimpleVector<Scalar>::extract(Scalar *v) const
    {
      memcpy(v, this->v, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->get_size() == vec->get_size());
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec->get(i));
      return this;

    }

    template<typename Scalar>
    Vector<Scalar>* SimpleVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->get_size(); i++)
        this->add(i, vec[i]);
      return this;
    }

    template<>
    HERMES_API SparseMatrix<double>* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new CSCMatrix<double>;
        }

      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new CSCMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new CSCMatrix<double>;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return NULL;
    }

    template<>
    HERMES_API Vector<double>* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new SimpleVector<double>;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<double>;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<double>;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionVector<double>;
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SimpleVector<double>;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return NULL;
    }

    template<>
    HERMES_API SparseMatrix<std::complex<double> >* create_matrix(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new CSCMatrix<std::complex<double> >;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new MumpsMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new CSCMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new CSCMatrix<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_matrix().");
      }
      return NULL;
    }

    template<>
    HERMES_API Vector<std::complex<double> >* create_vector(bool use_direct_solver)
    {
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          return new SimpleVector<std::complex<double> >;
        }

      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");

#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          return new EpetraVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          return new PetscVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION works only for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          return new SimpleVector<std::complex<double> >;
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_vector().");
      }
      return NULL;
    }

    template class Vector<double>;
    template class Vector<std::complex<double> >;

    template class SimpleVector<double>;
    template class SimpleVector<std::complex<double> >;
  }
}
