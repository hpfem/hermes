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
/*! \file superlu_solver.h
\brief SuperLU solver interface.
*/
#ifndef __HERMES_COMMON_SUPERLU_SOLVER_H_
#define __HERMES_COMMON_SUPERLU_SOLVER_H_

#include "config.h"
#ifdef WITH_SUPERLU
typedef int int_t;
#include "linear_matrix_solver.h"
#include "matrix.h"

#include <supermatrix.h>
#include <slu_util.h>
namespace Hermes
{
  namespace Solvers
  {
    template <typename Scalar> class SuperLUSolver;
#ifdef SLU_MT
    template <typename Scalar>
    class SuperLu
    {
    public:
      void gsequ (SuperMatrix *A, double *r, double *c, double *rowcnd, double *colcnd, double *amax, int *info);
      void laqgs (SuperMatrix *A, float *r, float *c, float rowcnd, float colcnd, float amax, char *equed);
      int gstrf (superlu_options_t *options, int m, int n, double anorm, LUstruct_t *LUstruct, gridinfo_t *grid, SuperLUStat_t *stat, int *info);
      float pivotGrowth (int ncols, SuperMatrix *A, int *perm_c, SuperMatrix *L, SuperMatrix *U);
      float langs (char *norm, SuperMatrix *A);
      void  gscon (char *norm, SuperMatrix *L, SuperMatrix *U, float anorm, float *rcond, SuperLUStat_t *stat, int *info);
      void  gstrs (trans_t trans, SuperMatrix *L, SuperMatrix *U, int *perm_c, int *perm_r, SuperMatrix *B, SuperLUStat_t *stat, int *info);
      double lamch_ (char *cmach);
      int querySpace (SuperMatrix *, SuperMatrix *, slu_memusage_t *);
    }
#else //SLU_MT

    typedef superlu_options_t         slu_options_t;
    typedef SuperLUStat_t             slu_stat_t;
    typedef struct
    {
      float for_lu;
      float total_needed;
    } slu_memusage_t;
#define SLU_DESTROY_L             Destroy_SuperNode_Matrix
#define SLU_DESTROY_U             Destroy_CompCol_Matrix
#define SLU_INIT_STAT(stat_ptr)   StatInit(stat_ptr)
#define SLU_PRINT_STAT(stat_ptr)  StatPrint(stat_ptr)

#define SLU_DTYPE                 SLU_D

#define SLU_PRINT_CSC_MATRIX    zPrint_CompCol_Matrix
#define Scalar_MALLOC       doublecomplexMalloc

    /** Type for storing number in SuperLU structures */
    template<typename Scalar> struct SuperLuType;

    /** Type for storing number in SuperLU real structures */
    template<>
    struct SuperLuType<double>
    {
    /** Type for storing scalar number in SuperLU real structures */
      typedef double Scalar;
    };

    /** Type for storing number in SuperLU complex structures */
    template<>
    struct SuperLuType<std::complex<double> >
    {
    /** Type for storing scalar number in SuperLU complex structures */
      typedef struct { double r, i; } Scalar;
    };
#endif //SLU_MT
  }
}

namespace Hermes
{
  namespace Algebra
  {
    /** \brief Matrix used with SuperLU solver */
    template <typename Scalar>
    class SuperLUMatrix : public SparseMatrix<Scalar>
    {
    public:
      SuperLUMatrix();
      virtual ~SuperLUMatrix();

      virtual void alloc();
      virtual void free();
      virtual Scalar get(unsigned int m, unsigned int n);
      virtual void zero();
      virtual void add(unsigned int m, unsigned int n, Scalar v);
      virtual void add_to_diagonal(Scalar v);
      virtual void add(unsigned int m, unsigned int n, Scalar **mat, int *rows, int *cols);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);
      virtual unsigned int get_matrix_size() const;
      virtual unsigned int get_nnz() const;
      virtual double get_fill_in() const;
      /// Add matrix.
      /// @param[in] mat matrix to be added
      virtual void add_matrix(SuperLUMatrix* mat);
      /// Add matrix to diagonal.
      /// @param[in] num_stages matrix is added to num_stages positions. num_stages * size(added matrix) = size(target matrix)
      /// @param[in] mat added matrix
      virtual void add_to_diagonal_blocks(int num_stages, SuperLUMatrix* mat);
      virtual void add_sparse_to_diagonal_blocks(int num_stages, SparseMatrix<Scalar>* mat);
      /// Add matrix to specific position.
      /// @param[in] i row in target matrix coresponding with top row of added matrix
      /// @param[in] j column in target matrix coresponding with lef column of added matrix
      /// @param[in] mat added matrix
      virtual void add_as_block(unsigned int i, unsigned int j, SuperLUMatrix* mat);

      // Applies the matrix to vector_in and saves result to vector_out.
      void multiply_with_vector(Scalar* vector_in, Scalar* vector_out);
      // Multiplies matrix with a Scalar.
      void multiply_with_Scalar(Scalar value);
      /// Creates matrix in SuperLU format using size, nnz, and the three arrays.
      /// @param[in] size size of matrix (num of rows and columns)
      /// @param[in] nnz number of nonzero values
      /// @param[in] ap index to ap/ax, where each column starts (size is matrix size + 1)
      /// @param[in] ai row indices
      /// @param[in] ax values
      void create(unsigned int size, unsigned int nnz, int* ap, int* ai, Scalar* ax);
      // Duplicates a matrix (including allocation).
      SuperLUMatrix<Scalar>* duplicate();

    protected:
      // SUPERLU specific data structures for storing the matrix (CSC format).
      /// Matrix entries (column-wise).
      Scalar *Ax;
      /// Row indices of values in Ax.
      int *Ai;
      /// Index to Ax/Ai, where each column starts.
      unsigned int *Ap;
      /// Number of non-zero entries ( =  Ap[size]).
      unsigned int nnz;

      friend class Solvers::SuperLUSolver<Scalar>;
      template<typename T> friend SparseMatrix<T>*  create_matrix();
    };

    /** \brief Vector used with SuperLU solver */
    template <typename Scalar>
    class SuperLUVector : public Vector<Scalar>
    {
    public:
      SuperLUVector();
      virtual ~SuperLUVector();

      virtual void alloc(unsigned int ndofs);
      virtual void free();
      virtual Scalar get(unsigned int idx);
      virtual void extract(Scalar *v) const;
      virtual void zero();
      virtual void change_sign();
      virtual void set(unsigned int idx, Scalar y);
      virtual void add(unsigned int idx, Scalar y);
      virtual void add(unsigned int n, unsigned int *idx, Scalar *y);
      virtual void add_vector(Vector<Scalar>* vec);
      virtual void add_vector(Scalar* vec);
      virtual bool dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt = DF_MATLAB_SPARSE);

    protected:
      /// SUPERLU specific data structures for storing the rhs.
      Scalar *v;     // Vector entries.

      friend class Solvers::SuperLUSolver<Scalar>;
    };
  }
  namespace Solvers
  {
    /// Encapsulation of SUPERLU linear solver.
    ///
    /// @ingroup solvers
    template <typename Scalar>
    class HERMES_API SuperLUSolver : public DirectSolver<Scalar>
    {
    public:
      /// Constructor of SuperLU solver.
      /// @param[in] m pointer to matrix
      /// @param[in] rhs pointer to right hand side vector
      SuperLUSolver(SuperLUMatrix<Scalar> *m, SuperLUVector<Scalar> *rhs);
      virtual ~SuperLUSolver();

      virtual bool solve();
      virtual int get_matrix_size();

    protected:
      /// Matrix to solve.
      SuperLUMatrix<Scalar> *m;
      /// Right hand side vector.
      SuperLUVector<Scalar> *rhs;

      bool has_A, has_B;            ///<  Have the native SuperLU matrices been created?
      bool inited;                  ///< Have the factorization structures been allocated?
      bool A_changed;               ///< Indicates that the system matrix has been changed.
      // internally during factorization or externally by
      // the user.

      /// Check the status returned from the solver routine.
      /// @param[in] info info returned from the solver routine
      bool check_status(unsigned int info);

      /// Deep copies of matrix and rhs data vectors (they may be changed by the solver driver,
      /// hence we need a copy so that the original SuperLUMatrix/Vector is preserved).
      int *local_Ai, *local_Ap;
      typename SuperLuType<Scalar>::Scalar *local_Ax, *local_rhs;

      bool setup_factorization();
      void free_factorization_data();
      void free_matrix();
      void free_rhs();

      SuperMatrix A, B;             ///< Native SuperLU representations of 'm' and 'rhs'.
      SuperMatrix L, U;             ///< L/U factors of A.
      double *R, *C;                ///< Row/column scaling factors of A.
      int *perm_r;                  ///< Row permutations from partial pivoting.
      int *perm_c;                  ///< Column permutations to reduce fill-in ( = > matrix Pc).
      int *etree;                   ///< Elimination tree of Pc'*A'*A*Pc.
      slu_options_t options;        ///< Structure holding the input options for the solver.

      private:
#ifndef SLU_MT
      void create_csc_matrix (SuperMatrix *A, int m, int n, int nnz, typename SuperLuType<Scalar>::Scalar *nzval, int *rowind, int *colptr,
        Stype_t stype, Dtype_t dtype, Mtype_t mtype);
      void  solver_driver (superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r, int *etree, char *equed, double *R,
        double *C, SuperMatrix *L, SuperMatrix *U, void *work, int lwork, SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth,
        double *rcond, double *ferr, double *berr, slu_memusage_t *mem_usage, SuperLUStat_t *stat, int *info);
      void create_dense_matrix (SuperMatrix *X, int m, int n, typename SuperLuType<Scalar>::Scalar *x, int ldx, Stype_t stype, Dtype_t dtype, Mtype_t mtype);
#endif  //SLU_MT

#ifndef SLU_MT
      char equed[1];              ///< Form of equilibration that was done on A.
#else
      equed_t equed;              ///< Form of equilibration that was done on A.
      SuperMatrix AC;             ///< Matrix A permuted by perm_c.
#endif //SLU_MT
      template<typename T> friend LinearMatrixSolver<T>* create_linear_solver(Matrix<T>* matrix, Vector<T>* rhs);
    };
  }
}
#endif //SUPER_LU
#endif