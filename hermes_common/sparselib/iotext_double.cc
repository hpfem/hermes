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

#include<stdio.h>
#include<stdlib.h>
#include <iostream>

#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"

void readtxtfile_mat(const char *filename, Coord_Mat_double *A)
{
/*************************************************************************/
/*  This function opens and reads the specified file, interpreting its   */
/*  contents as a sparse matrix stored one nonzero per line in the       */
/*  format:                                                              */
/*                                                                       */
/*  i   j    A(i,j)                                                      */
/*                                                                       */
/*  Although the internal storage is 0-index based, the I/O functions    */  
/*  are 1-index based for general Fortran compatibility.                 */  
/*                                                                       */
/*    ----------                                                         */
/*    **CAVEAT**                                                         */
/*    ----------                                                         */
/*  To determine the appropriate matrix size, this function assumes that */
/*  the A(A.dim(0),A.dim(1)) element appears in the file, even if it is  */
/*  a zero element.                                                      */
/*                                                                       */
/*************************************************************************/
    FILE *in_file;
    char line[82];
    char* line_ptr;
    
  
    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       std::cerr << "Cannot open file: " << filename << "\n";
       exit(1);
    }
    
//  Read through file first, counting nonzero elements and looking for
//  matrix dimensions...

    int args, i, j;
    double value;
    int count = 0;
    int maxrow = 0;
    int maxcol = 0;
    while (1)
    {
       line_ptr = fgets(line, 82, in_file);
       if (line_ptr == NULL) break;
       args = sscanf(line_ptr,"%d %d %le",&i,&j,&value);
       if (args != 3) 
       {
          printf("Error reading textfile:%s\n",filename);
          exit(1);
       }
       if (i > maxrow) maxrow = i;
       if (j > maxcol) maxcol = j;
       count++;
    }
    fclose(in_file);

//  create the arrays to hold the file information:
    
    double *val= new double[count];
    int *colind = new int[count], *rowind = new int[count];

    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       std::cerr << "Cannot open file: " << filename << "\n";
       exit(1);
    }
    
    for (i=0;i<count;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       if (line_ptr == NULL) break;
       args = sscanf(line_ptr,"%d %d %le",&rowind[i],&colind[i],&val[i]);
       rowind[i]--; colind[i]--;
       if (args != 3)
       {
          printf("Error reading textfile:%s\n",filename);
          exit(1);
       }
    }

    Coord_Mat_double C(maxrow,maxcol,count, val, rowind, colind);
    *A = C;

    return;
}


void readtxtfile_mat(const char *filename, CompCol_Mat_double *A)
{
    Coord_Mat_double C;
    readtxtfile_mat(filename, &C);
    *A = C;
    return;
}

void readtxtfile_mat(const char *filename, CompRow_Mat_double *A)
{
    Coord_Mat_double C;
    readtxtfile_mat(filename, &C);
    *A = C;
    return;
}

void writetxtfile_mat(const char *filename, const Coord_Mat_double &A)
{
/*************************************************************************/
/*  This function opens and writes to specified file the nonzero entries */
/*  of A.   To ensure that the row and column dimensions, M and N, of    */
/*  the matrix can be determined implicitly from the file, we include in */
/*  the output the A(M-1,N-1) matrix element, even if it falls outside   */
/*  of the sparsity pattern of A.  Each row of the output file has the   */
/*  format:                                                              */
/*                                                                       */
/*  i   j    A(i,j)                                                      */
/*                                                                       */
/*  Although the internal storage is 0-index based, the I/O functions    */
/*  are 1-index based for general Fortran compatibility.                 */
/*                                                                       */
/*************************************************************************/
    FILE *out_file;
    out_file = fopen( filename, "w");

    int nnz = A.NumNonzeros();
    int rowp1, colp1;
    int M = A.dim(0);
    int N = A.dim(1);
    int flag = 0;
//  Loop through Nonzeros
    for (int j = 0; j < nnz ; j++)
    {
       rowp1 = A.row_ind(j) +1;
       colp1 = A.col_ind(j) +1;
       if ( rowp1 == M && colp1 == N ) flag = 1;
       fprintf(out_file, "%14d\t%14d\t%20.16e\n", rowp1, colp1,A.val(j));
    }
    if (flag == 0)
       fprintf(out_file, "%14d\t%14d\t%20.16e\n", M, N, A(M-1,N-1));
    fclose(out_file);

    return;
}


void writetxtfile_mat(const char *filename, const CompCol_Mat_double &A)
{
    FILE *out_file;
    out_file = fopen( filename, "w");
 
    int rowp1, colp1;
    int M = A.dim(0);
    int N = A.dim(1);
    int flag = 0;
//  Loop through columns
    for (int j = 0; j < N ; j++)
       for (int i=A.col_ptr(j);i<A.col_ptr(j+1);i++)
       {
          rowp1 = A.row_ind(i)+1;
          colp1 = j + 1;
          if ( rowp1 == M && colp1 == N ) flag = 1;
          fprintf(out_file,"%14d%4s%14d%4s% 20.16E\n", rowp1, "    ",
                                                       colp1,"    ", A.val(i));
       }
 
    if (flag == 0)
       fprintf(out_file,"%14d\t%14d\t%20.16E\n", M,  N, A(M-1,N-1));

    fclose(out_file);

    return;
}

void writetxtfile_mat(const char *filename, const CompRow_Mat_double &A)
{
    FILE *out_file;
    out_file = fopen( filename, "w");
 
    int rowp1, colp1;
    int M = A.dim(0);
    int N = A.dim(1);
    int flag = 0;
//  Loop through rows...
    for (int i = 0; i < M ; i++)
       for (int j=A.row_ptr(i);j<A.row_ptr(i+1);j++)
       {
          rowp1 =  i + 1;
          colp1 =  A.col_ind(j) + 1;
          if ( rowp1 == M && colp1 == N ) flag = 1;
          fprintf(out_file,"%14d\t%14d\t%20.16e\n", rowp1, colp1,A.val(j));
       }
 
    if (flag == 0)
       fprintf(out_file,"%14d\t%14d\t%20.16e\n", M, N, A(M-1,N-1));

    fclose(out_file);

    return;
}

