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


#ifndef _IOHB_H_
#define _IOHB_H_

#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif
 
#ifdef __cplusplus
extern "C" {
#endif

void  readHB_info       __ProtoGlarp__(( const char*, int*, int*, int*, int* ));
void readHB_mat_double  __ProtoGlarp__(( const char*, int*, int*, double*));
void readHB_mat_float   __ProtoGlarp__(( const char*, int*, int*, float*));

void readHB_rhs_double __ProtoGlarp__(( const char*, double*, int));
void readHB_rhs_float  __ProtoGlarp__(( const char*, float*, int));

void writeHB_mat_double 
    __ProtoGlarp__(( const char*, int, int, int, const int*, const int*, 
            const double*, int, const double*, const char*, const char*));
void writeHB_mat_float  
    __ProtoGlarp__(( const char*, int, int, int, const int*, const int*, 
            const float*, int, const double*, const char*, const char*));

#ifdef __cplusplus
}
#endif


#endif
