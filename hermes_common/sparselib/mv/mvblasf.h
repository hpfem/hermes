
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

#ifndef _MV_BLAS1_float_H_
#define _MV_BLAS1_float_H_


MV_Vector_float& operator*=(MV_Vector_float &x, const float &a);
MV_Vector_float operator*(const float &a, const MV_Vector_float &x);
MV_Vector_float operator*(const MV_Vector_float &x, const float &a);
MV_Vector_float operator+(const MV_Vector_float &x, 
    const MV_Vector_float &y);
MV_Vector_float operator-(const MV_Vector_float &x, 
    const MV_Vector_float &y);
MV_Vector_float& operator+=(MV_Vector_float &x, const MV_Vector_float &y);
MV_Vector_float& operator-=(MV_Vector_float &x, const MV_Vector_float &y);

float dot(const MV_Vector_float &x, const MV_Vector_float &y);
float norm(const MV_Vector_float &x);

#endif

// _MV_BLAS1_float_H_
