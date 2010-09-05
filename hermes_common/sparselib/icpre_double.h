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

#ifndef ICPRE_H
#define ICPRE_H

#include "vecdefs.h"
#include VECTOR_H
#include "comprow_double.h"
#include "compcol_double.h"


class ICPreconditioner_double {

private:
  VECTOR_double val_;
  VECTOR_int    pntr_;
  VECTOR_int    indx_;
  int nz_;
  int dim_[2];

public:
  ICPreconditioner_double(const CompCol_Mat_double &A);
  ICPreconditioner_double(const CompRow_Mat_double &A);
  ~ICPreconditioner_double(void){};

  VECTOR_double     solve(const VECTOR_double &x) const;
  VECTOR_double     trans_solve(const VECTOR_double &x) const;
};

#endif
