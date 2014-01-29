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
/*! \file mumps_solver.h
\brief MUMPS solver interface.
*/
#ifndef __HERMES_COMMON_MUMPS_SOLVER_H_
#define __HERMES_COMMON_MUMPS_SOLVER_H_
#include "config.h"
#ifdef WITH_MUMPS
#include "solvers/linear_matrix_solver.h"
#include "algebra/matrix.h"

extern "C"
{
#include <mumps_c_types.h>
#include <dmumps_c.h>
#include <zmumps_c.h>
}

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class MumpsSolver;
  }
}
namespace Hermes
{
  namespace Algebra
  {
    /** Type for storing number in Mumps structures */
    template <typename Scalar> struct mumps_type;

    /** Type for storing number in Mumps complex structures */
    template <>
    struct mumps_type<std::complex<double> >
    {
      /** Type for storing mumps struct in Mumps complex structures */
      typedef ZMUMPS_STRUC_C mumps_struct;
      /** Type for storing scalar number in Mumps complex structures */
      typedef ZMUMPS_COMPLEX mumps_Scalar;
    };

    /** Type for storing number in Mumps real structures */
    template <>
    struct mumps_type<double>
    {
      /** Type for storing mumps struct in Mumps real structures */
      typedef DMUMPS_STRUC_C mumps_struct;
      /** Type for storing scalar number in Mumps real structures */
      typedef double mumps_Scalar;
    };

    /** \brief Matrix used with MUMPS solver */
    /// Important: MUMPS is indexing from 1
    template <typename Scalar>
    class MumpsMatrix : public CSCMatrix<Scalar>
    {
    public:
      MumpsMatrix();
      virtual ~MumpsMatrix();

      void alloc_data();
      void free();
      Scalar get(unsigned int m, unsigned int n) const;
      void zero();

      void add(unsigned int m, unsigned int n, Scalar v);

      /// Matrix export method.
      /// Utility version
      /// \See MatrixRhsImportExport<Scalar>::export_to_file.
      void export_to_file(const char *filename, const char *var_name, MatrixExportFormat fmt, char* number_format = "%lf");

      /// Reading matrix
      /// Utility version
      /// \See Matrix<Scalar>::import_from_file.
      void import_from_file(const char *filename, const char *var_name, MatrixExportFormat fmt);

      /// Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);

      /// Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar*& vector_out, bool vector_out_initialized = false) const;

      /// Creates matrix using size, nnz, and the three arrays.
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);

      /// Duplicates a matrix (including allocation).
      CSMatrix<Scalar>* duplicate() const;

    protected:
      int *irn;         ///< Row indices.
      int *jcn;         ///< Column indices.
      typename mumps_type<Scalar>::mumps_Scalar *Ax; ///< Matrix entries (column-wise).

      friend class Solvers::MumpsSolver<Scalar>;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };
  }

  namespace Solvers
  {
    /// Encapsulation of MUMPS linear solver.
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API MumpsSolver : public DirectSolver<Scalar>
    {
    public:
      /// Constructor of MumpsSolver.
      /// @param[in] m matrix pointer
      /// @param[in] rhs right hand side pointer
      MumpsSolver(MumpsMatrix<Scalar> *m, SimpleVector<Scalar> *rhs);
      virtual ~MumpsSolver();
      void free();

      virtual void solve();
      virtual int get_matrix_size();

      /// Matrix to solve.
      MumpsMatrix<Scalar> *m;
      /// Right hand side.
      SimpleVector<Scalar> *rhs;

      /// \todo document
      bool setup_factorization();

      /// MUMPS specific structure with solver parameters.
      typename mumps_type<Scalar>::mumps_struct param;

      /// @return false on error
      bool check_status();

      /// (Re)initialize a MUMPS instance.
      /// @return true on succes
      /// \sa #check_status()
      bool reinit();
    private:
      void mumps_c(typename mumps_type<Scalar>::mumps_struct * param);  //wrapper around dmums_c or zmumps_c

      /// True if solver is inited.
      bool inited;

      /// Internal - control parameter for MUMPS.
      /// See MUMPS doc, page 27, version 4.10.
      int icntl_14;
      /// Initial value is 1 * 100%.
      static const int init_icntl_14 = 1;
      /// Maximum value is 128 * 100%.
      static const int max_icntl_14 = 128;
    };
  }
}
#endif
#endif
