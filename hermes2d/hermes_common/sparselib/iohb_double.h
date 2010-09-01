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
/*      Read a Harwell-Boeing file into a compressed Column matrix       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#include "iohb.h"


CompCol_Mat_double& readHB_mat(const char *filename, CompCol_Mat_double *A);
CompRow_Mat_double& readHB_mat(const char *filename, CompRow_Mat_double *A);
Coord_Mat_double& readHB_mat(const char *filename, Coord_Mat_double *A);

const CompCol_Mat_double& writeHB_mat(const char *filename, 
    const CompCol_Mat_double &A, 
    int nrhs=0, const double *rhs = 0, 
    const char *title=0, const char *key=0);

const CompRow_Mat_double& writeHB_mat(const char *filename, 
    const CompRow_Mat_double &A, 
    int nrhs=0, const double* rhs = 0,
    const char *title=0, const char *key=0);


const Coord_Mat_double& writeHB_mat(const char *filename, 
    const Coord_Mat_double &A, 
    int nrhs=0, const double* rhs=0,
    const char *title=0, const char *key=0);

VECTOR_double& readHB_rhs(const char *filename, VECTOR_double *b, int j=0);
