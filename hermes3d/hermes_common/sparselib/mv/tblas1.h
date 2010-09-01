
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

#ifndef _MV_BLAS1_TYPE_H_
#define _MV_BLAS1_TYPE_H_

#include <math.h>
#include <stdlib.h>


MV_Vector_TYPE& operator*=(MV_Vector_TYPE &x, const TYPE &a);
MV_Vector_TYPE operator*(const TYPE &a, const MV_Vector_TYPE &x);
MV_Vector_TYPE operator*(const MV_Vector_TYPE &x, const TYPE &a);
MV_Vector_TYPE operator+(const MV_Vector_TYPE &x, 
    const MV_Vector_TYPE &y);
MV_Vector_TYPE operator-(const MV_Vector_TYPE &x, 
    const MV_Vector_TYPE &y);
MV_Vector_TYPE& operator+=(MV_Vector_TYPE &x, const MV_Vector_TYPE &y);
MV_Vector_TYPE& operator-=(MV_Vector_TYPE &x, const MV_Vector_TYPE &y);

TYPE dot(const MV_Vector_TYPE &x, const MV_Vector_TYPE &y);
TYPE norm(const MV_Vector_TYPE &x);

#endif
// _MV_BLAS1_TYPE_H_
