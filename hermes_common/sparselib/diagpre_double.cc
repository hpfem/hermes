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

#include <iostream>
#include <stdlib.h>
#include "diagpre_double.h"

int 
CopyInvDiagonals(int n, const int *pntr, const int *indx,
         const double *sa, double *diag)
{
  int i, j;

  for (i = 0; i < n; i++)
    diag[i] = 0;
  
  /* Find the diagonal elements */
  for (i = 0; i < n; i++) {
    for (j = pntr[i]; j < pntr[i+1]; j++) {
      if (indx[j] == i) {
    if (sa[j] == 0)
      return i;
    diag[i] = 1. / sa[j];
    break;
      }
    }
    if (diag[i] == 0)
      return -i;
  }
  
  return 0;
}


DiagPreconditioner_double::DiagPreconditioner_double(const CompRow_Mat_double &C)
  : diag_(C.dim(0))
{
  int i;

  if ((i = CopyInvDiagonals(C.dim(0), &C.row_ptr(0), &C.col_ind(0), 
                &C.val(0), &diag(0))) != 0) {
    std::cerr << "Diagonal preconditioner failure.";
    std::cerr << " Zero detected in element " << i << "\n";
    exit(1);
  }
}


DiagPreconditioner_double::DiagPreconditioner_double(const CompCol_Mat_double &C)
  : diag_ (C.dim (0))
{
  int i;

  if ((i = CopyInvDiagonals(C.dim(0), &C.col_ptr(0), &C.row_ind(0), 
                &C.val(0), &diag(0))) != 0) {
    std::cerr << "Diagonal preconditioner failure.";
    std::cerr << " Zero detected in element " << i << "\n";
    exit(1);
  }
}


VECTOR_double 
DiagPreconditioner_double::solve (const VECTOR_double &x) const 
{
  VECTOR_double y(x.size());

  for (int i = 0; i < x.size(); i++)
    y(i) = x(i) * diag(i);
  
  return y;
}


VECTOR_double
DiagPreconditioner_double::trans_solve (const VECTOR_double &x) const 
{
  VECTOR_double y(x.size());

  for (int i = 0; i < x.size(); i++)
    y(i) = x(i) * diag(i);
  
  return y;
}
