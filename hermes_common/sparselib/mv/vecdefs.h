
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*      Which dense vector/matrix classes to build SparseLib++ from      */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef vector_defs_H
#define vector_defs_H

#define VECTOR_H              "mvv.h"
#define VECTOR_double         MV_Vector_double
#define VECTOR_float          MV_Vector_float
#define VECTOR_int            MV_Vector_int
#define VECTOR_ref            MV_Vector_::ref
#define MATRIX_H              "mvm.h"
#define MATRIX_double         MV_ColMat_double
#define MATRIX_float          MV_ColMat_float
#define MATRIX_int            MV_ColMat_int
#define MATRIX_ref            MV_Matrix_::ref

#define VECTOR_COMPLEX        MV_Vector_COMPLEX
#define MATRIX_COMPLEX        MV_ColMat_COMPLEX

#endif
