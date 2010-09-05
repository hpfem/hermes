
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

#ifndef _MV_BLAS1_double_H_
#define _MV_BLAS1_double_H_


MV_Vector_double& operator*=(MV_Vector_double &x, const double &a);
MV_Vector_double operator*(const double &a, const MV_Vector_double &x);
MV_Vector_double operator*(const MV_Vector_double &x, const double &a);
MV_Vector_double operator+(const MV_Vector_double &x, 
    const MV_Vector_double &y);
MV_Vector_double operator-(const MV_Vector_double &x, 
    const MV_Vector_double &y);
MV_Vector_double& operator+=(MV_Vector_double &x, const MV_Vector_double &y);
MV_Vector_double& operator-=(MV_Vector_double &x, const MV_Vector_double &y);

double dot(const MV_Vector_double &x, const MV_Vector_double &y);
double norm(const MV_Vector_double &x);

#endif

// _MV_BLAS1_double_H_
