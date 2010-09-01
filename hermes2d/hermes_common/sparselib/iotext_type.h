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

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*      I/O for plain ascii text files                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "compcol_TYPE.h"
#include "comprow_TYPE.h"
#include "coord_TYPE.h"

void readtxtfile_mat(const char *filename, Coord_Mat_TYPE &A);
void readtxtfile_mat(const char *filename, CompCol_Mat_TYPE &A);
void readtxtfile_mat(const char *filename, CompRow_Mat_TYPE &A);
void writetxtfile_mat(const char *filename, Coord_Mat_TYPE &A);
void writetxtfile_mat(const char *filename, CompCol_Mat_TYPE &A);
void writetxtfile_mat(const char *filename, CompRow_Mat_TYPE &A);
