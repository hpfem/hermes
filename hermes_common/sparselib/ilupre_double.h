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

#ifndef ILUPRE_H
#define ILUPRE_H

#include "vecdefs.h"
#include VECTOR_H
#include "comprow_double.h"
#include "compcol_double.h"


class CompCol_ILUPreconditioner_double {

 private:
  VECTOR_double l_val_;
  VECTOR_int    l_colptr_;
  VECTOR_int    l_rowind_;
  int l_nz_;

  VECTOR_double u_val_;
  VECTOR_int    u_colptr_;
  VECTOR_int    u_rowind_;
  int u_nz_;

  int dim_[2];
  
 public:
  CompCol_ILUPreconditioner_double(const CompCol_Mat_double &A);
  ~CompCol_ILUPreconditioner_double(void){};
  
  VECTOR_double     solve(const VECTOR_double &x) const;
  VECTOR_double     trans_solve(const VECTOR_double &x) const;
};


class CompRow_ILUPreconditioner_double {

 private:
  VECTOR_double l_val_;
  VECTOR_int    l_rowptr_;
  VECTOR_int    l_colind_;
  int l_nz_;

  VECTOR_double u_val_;
  VECTOR_int    u_rowptr_;
  VECTOR_int    u_colind_;
  int u_nz_;

  int dim_[2];

 public:
  CompRow_ILUPreconditioner_double(const CompRow_Mat_double &A);
  ~CompRow_ILUPreconditioner_double(void){};
  
  VECTOR_double     solve(const VECTOR_double &x) const;
  VECTOR_double     trans_solve(const VECTOR_double &x) const;
};


#endif
