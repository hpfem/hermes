/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***                              */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdlib.h>
#include "ilupre_double.h"
#include "qsort_double.h"

#include "spblas.h"

CompCol_ILUPreconditioner_double::
CompCol_ILUPreconditioner_double(const CompCol_Mat_double &A)
: 
  l_val_(0), l_colptr_(A.dim(1) + 1), l_rowind_(0), l_nz_(0), u_val_(0),
    u_colptr_(A.dim(1) + 1), u_rowind_(0), u_nz_(0)
{
  int i, j, k, pn, qn, rn;
  double multiplier;

  // Copy
  dim_[0] = A.dim(0);
  dim_[1] = A.dim(1);


  // Get size of l and u
  for (i = 0; i < dim_[1]; i++)
    for (j = A.col_ptr(i); j < A.col_ptr(i+1); j++)
      if (A.row_ind(j) > i)
    l_nz_++;
      else
    u_nz_++;

  l_val_.newsize(l_nz_);
  u_val_.newsize(u_nz_);
  l_rowind_.newsize(l_nz_);
  u_rowind_.newsize(u_nz_);

  l_colptr_(0) = u_colptr_(0) = 0;

  // Split up A into l and u
  for (i = 0; i < dim_[1]; i++) {
    l_colptr_(i+1) = l_colptr_(i);
    u_colptr_(i+1) = u_colptr_(i);
    
    for (j = A.col_ptr(i); j < A.col_ptr(i+1); j++)
      if (A.row_ind(j) > i) {
    k = l_colptr_(i+1)++;
    l_val_(k) = A.val(j);
    l_rowind_(k) = A.row_ind(j);
      } else if (A.row_ind(j) <= i) {
    k = u_colptr_(i+1)++;
    u_val_(k) = A.val(j);
    u_rowind_(k) = A.row_ind(j);
      }
  }
  
  for (i = 0; i < dim_[1]; i++) {
    QSort(l_rowind_, l_val_, l_colptr_[i], l_colptr_[i+1] - l_colptr_[i]);
    QSort(u_rowind_, u_val_, u_colptr_[i], u_colptr_[i+1] - u_colptr_[i]);
  }

  // Factor matrix
  for (i = 0; i < dim_[0] - 1; i++) {
    multiplier = u_val_(u_colptr_(i+1)-1);
    
    for (j = l_colptr_[i]; j < l_colptr_[i+1]; j++)
      l_val_[j] /= multiplier;
    
    for (j = u_colptr_[i+1]; j < u_colptr_[i+2]-1; j++) {
      multiplier = u_val_[j];
      qn = j + 1;
      rn = l_colptr_[i+1];
      for (pn = l_colptr_[u_rowind_[j]]; 
       pn < l_colptr_[u_rowind_[j]+1] && l_rowind_[pn] <= i + 1; 
       pn++) {
    while (qn < u_colptr_[i+2] && u_rowind_[qn] < l_rowind_[pn])
      qn++;
    if (qn < u_colptr_[i+2] && l_rowind_[pn] == u_rowind_[qn])
      u_val_[qn] -= multiplier * l_val_[pn];
      }
      for ( ; pn < l_colptr_[u_rowind_[j]+1]; pn++) {
    while (rn < l_colptr_[i+2] && l_rowind_[rn] < l_rowind_[pn])
      rn++;
    if (rn < l_colptr_[i+2] && l_rowind_[pn] == l_rowind_[rn])
      l_val_[rn] -= multiplier * l_val_[pn];
      }
    }
  }
}


CompRow_ILUPreconditioner_double::
CompRow_ILUPreconditioner_double(const CompRow_Mat_double &A)
: l_val_(0), l_rowptr_(A.dim(1) + 1), l_colind_(0), l_nz_(0),
    u_val_(0), u_rowptr_(A.dim(1) + 1), u_colind_(0), u_nz_(0)
{
  int i, j, k, pn, qn, rn;
  double multiplier;

  // Copy
  dim_[0] = A.dim(0);
  dim_[1] = A.dim(1);


  // Get size of l and u
  for (i = 0; i < dim_[1]; i++)
    for (j = A.row_ptr(i); j < A.row_ptr(i+1); j++)
      if (A.col_ind(j) < i)
    l_nz_++;
      else
    u_nz_++;

  l_val_.newsize(l_nz_);
  u_val_.newsize(u_nz_);
  l_colind_.newsize(l_nz_);
  u_colind_.newsize(u_nz_);

  l_rowptr_(0) = u_rowptr_(0) = 0;

  // Split up A into l and u
  for (i = 0; i < dim_[1]; i++) {
    l_rowptr_(i+1) = l_rowptr_(i);
    u_rowptr_(i+1) = u_rowptr_(i);
    
    for (j = A.row_ptr(i); j < A.row_ptr(i+1); j++)
      if (A.col_ind(j) < i) {
    k = l_rowptr_(i+1)++;
    l_val_(k) = A.val(j);
    l_colind_(k) = A.col_ind(j);
      } else if (A.col_ind(j) >= i) {
    k = u_rowptr_(i+1)++;
    u_val_(k) = A.val(j);
    u_colind_(k) = A.col_ind(j);
      }
  }

  for (i = 0; i < dim_[1]; i++) {
    QSort(l_colind_, l_val_, l_rowptr_[i], l_rowptr_[i+1] - l_rowptr_[i]);
    QSort(u_colind_, u_val_, u_rowptr_[i], u_rowptr_[i+1] - u_rowptr_[i]);
  }

  // Factor matrix
  for (i = 1; i < dim_[0]; i++) {
    for (j = l_rowptr_[i]; j < l_rowptr_[i+1]; j++) {
      pn = u_rowptr_[l_colind_[j]];
      multiplier = (l_val_[j] /= u_val_[pn]);

      qn = j + 1;
      rn = u_rowptr_[i];

      for (pn++; pn < u_rowptr_[l_colind_[j]+1] && u_colind_[pn] < i; pn++) {
    while (qn < l_rowptr_[i+1] && l_colind_[qn] < u_colind_[pn])
      qn++;
    if (qn < l_rowptr_[i+1] && u_colind_[pn] == l_colind_[qn])
      l_val_[qn] -= multiplier * u_val_[pn];
      }
      for ( ; pn < u_rowptr_[l_colind_[j]+1]; pn++) {
    while (rn < u_rowptr_[i+1] && u_colind_[rn] < u_colind_[pn])
      rn++;
    if (rn < u_rowptr_[i+1] && u_colind_[pn] == u_colind_[rn])
      u_val_[rn] -= multiplier * u_val_[pn];
      }
    }
  }
}


VECTOR_double
CompCol_ILUPreconditioner_double::solve(const VECTOR_double &x) const 
{
  int M = x.size();
  VECTOR_double y(M);
  VECTOR_double work(M);

  int descra[9];

  descra[0] = 0;

  // lower unit
  descra[1] = 1;
  descra[2] = 1;

  F77NAME(dcscsm) (0, M, 1, 1, NULL, 1.0,
           descra, &l_val_(0), &l_rowind_(0), &l_colptr_(0),
           &x(0), M, 0.0, &y(1), M,
           &work(0), M);

  // upper diag
  descra[1] = 2;
  descra[2] = 0;

  F77NAME(dcscsm) (0, M, 1, 1, NULL, 1.0,
           descra, &u_val_(0), &u_rowind_(0), &u_colptr_(0),
           &y(0), M, 0.0, &y(1), M,
           &work(0), M);

  return y;
} 


VECTOR_double
CompCol_ILUPreconditioner_double::trans_solve(const VECTOR_double &x) const 
{
  int M = x.size();
  VECTOR_double y(M);
  VECTOR_double work(M);

  int descra[9];

  descra[0] = 0;

  // upper diag transpose
  descra[1] = 2;
  descra[2] = 0;

  F77NAME(dcscsm) (1, M, 1, 1, NULL, 1.0,
           descra, &u_val_(0), &u_rowind_(0), &u_colptr_(0),
           &x(0), M, 0.0, &y(1), M,
           &work(0), M);

  // lower unit transpose
  descra[1] = 1;
  descra[2] = 1;

  F77NAME(dcscsm) (1, M, 1, 1, NULL, 1.0,
           descra, &l_val_(0), &l_rowind_(0), &l_colptr_(0),
           &y(0), M, 0.0, &y(1), M,
           &work(0), M);

  return y;
}


VECTOR_double
CompRow_ILUPreconditioner_double::solve(const VECTOR_double &x) const 
{
  int M = x.size();
  VECTOR_double y(M);
  VECTOR_double work(M);

  int descra[9];

  descra[0] = 0;

  // lower unit
  descra[1] = 1;
  descra[2] = 1;

  F77NAME(dcsrsm) (0, M, 1, 1, NULL, 1.0,
           descra, &l_val_(0), &l_colind_(0), &l_rowptr_(0),
           &x(0), M, 0.0, &y(1), M,
           &work(0), M);

  // upper diag
  descra[1] = 2;
  descra[2] = 0;

  F77NAME(dcsrsm) (0, M, 1, 1, NULL, 1.0,
           descra, &u_val_(0), &u_colind_(0), &u_rowptr_(0),
           &y(0), M, 0.0, &y(1), M,
           &work(0), M);

  return y;
} 


VECTOR_double
CompRow_ILUPreconditioner_double::trans_solve(const VECTOR_double &x) const 
{
  int M = x.size();
  VECTOR_double y(M);
  VECTOR_double work(M);

  int descra[9];

  descra[0] = 0;

  // upper diag transpose
  descra[1] = 2;
  descra[2] = 0;

  F77NAME(dcsrsm) (1, M, 1, 1, NULL, 1.0,
           descra, &u_val_(0), &u_colind_(0), &u_rowptr_(0),
           &x(0), M, 0.0, &y(1), M,
           &work(0), M);

  // lower unit transpose
  descra[1] = 1;
  descra[2] = 1;

  F77NAME(dcsrsm) (1, M, 1, 1, NULL, 1.0,
           descra, &l_val_(0), &l_colind_(0), &l_rowptr_(0),
           &y(0), M, 0.0, &y(1), M,
           &work(0), M);

  return y;
} 
