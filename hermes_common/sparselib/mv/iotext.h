
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

//  read and write MV vectors as text files with (at most) one element
//  per line.

#ifndef _IOTEXT_H_
#define _IOTEXT_H_

int readtxtfile_vec(const char *filename, MV_Vector_double *Aptr);
int writetxtfile_vec(const char *filename, const MV_Vector_double &A);
int readtxtfile_vec(const char *filename, MV_Vector_int *Aptr);
int writetxtfile_vec(const char *filename, const MV_Vector_int &A);

#endif

