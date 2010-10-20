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
/*           Compressed row sparse matrix (0-based)                      */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <iostream>
#include <stdlib.h>
#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"

#include "spblas.h"

/*****************************/
/*  Constructor(s)           */
/*****************************/

CompRow_Mat_double::CompRow_Mat_double(void)
        : val_(0), rowptr_(0), colind_(0), base_(0), nz_(0)
{
        dim_[0] = 0;
        dim_[1] = 0;
}

/*****************************/
/*  Copy constructor         */
/*****************************/

CompRow_Mat_double::CompRow_Mat_double(const CompRow_Mat_double &S) :  
        val_(S.val_), rowptr_(S.rowptr_), colind_(S.colind_), base_(S.base_), 
        nz_(S.nz_)
{
        dim_[0] = S.dim_[0];
        dim_[1] = S.dim_[1];
}

/***********************************/
/*  Construct from storage vectors */
/***********************************/

CompRow_Mat_double::CompRow_Mat_double(int M, int N, int nz, double *val, 
                                     int *r, int *c, int base) :
        val_(val, nz), rowptr_(r, M+1), colind_(c, nz), base_(base), nz_(nz) 
{
        dim_[0] = M;
        dim_[1] = N;
}

CompRow_Mat_double::CompRow_Mat_double(int M, int N, int nz, 
                   const VECTOR_double &val, const VECTOR_int &r, 
                   const VECTOR_int &c, int base) :
        val_(val), rowptr_(r), colind_(c), base_(base), nz_(nz)
{
        dim_[0] = M;
        dim_[1] = N;
} 

/**********************************************************************/
/*  Construct a CompRow_Mat_double from a CompCol_Mat_double            */
/*  (i.e. convert compressed column storage to compressed row storage)*/
/**********************************************************************/

CompRow_Mat_double::CompRow_Mat_double(const CompCol_Mat_double &C) : 
         val_(C.NumNonzeros()), rowptr_(C.dim(0) +1), 
        colind_(C.NumNonzeros()), base_(C.base()), nz_(C.NumNonzeros())
{

        dim_[0] = C.dim(0);
        dim_[1] = C.dim(1);
        
        int i,j;

        VECTOR_int tally(C.dim(0)+1, 0);
//      First pass through nonzeros.  Tally entries in each row.
//      And calculate rowptr array.
        for (i=0;i<nz_;i++) tally(C.row_ind(i))++;
        rowptr_(0) = 0;
        for (j=0;j<dim_[0];j++) rowptr_(j+1) = rowptr_(j)+tally(j);
//      Make copy of rowptr for use in second pass.
        tally = rowptr_;
 
//      Second pass through nonzeros.   Fill in index and value entries.
        int count = 0;
        for (i=1;i<=dim_[1];i++)
        {
           for (j=count;j<C.col_ptr(i);j++)
           {
              val_(tally(C.row_ind(j))) = C.val(j);
              colind_(tally(C.row_ind(j))) = i-1;
              tally(C.row_ind(count))++;
              count++;
           }
        }
}

/**********************************************************************/
/*  Construct a CompRow_Mat_double from a Coord_Mat_double              */
/*  (i.e. convert coordinate storage to compressed row storage)       */
/**********************************************************************/

CompRow_Mat_double::CompRow_Mat_double(const Coord_Mat_double &CO) : 
        val_(CO.NumNonzeros()), rowptr_(CO.dim(0)+1), 
        colind_(CO.NumNonzeros()), base_(CO.base()), nz_(CO.NumNonzeros())
{

        dim_[0] = CO.dim(0);
        dim_[1] = CO.dim(1);

        int i;
        VECTOR_int tally(CO.dim(0)+1, 0);
//      First pass through nonzeros.  Tally entries in each row.
//      And calculate rowptr array.
        for (i=0;i<nz_;i++) tally(CO.row_ind(i))++;
        rowptr_(0) = 0;
        for (i=0;i<dim_[0];i++) rowptr_(i+1) = rowptr_(i)+tally(i);
//      Make copy of rowptr for use in second pass.
        tally = rowptr_;
//      Second pass through nonzeros.   Fill in index and value entries.
        for (i=0;i<nz_;i++)
        {
           val_(tally(CO.row_ind(i))) = CO.val(i);
           colind_(tally(CO.row_ind(i))) = CO.col_ind(i);
           tally(CO.row_ind(i))++;
        }
}
              

/***************************/
/* Assignment operator...  */
/***************************/

CompRow_Mat_double& CompRow_Mat_double::operator=(const CompRow_Mat_double &R)
{
        dim_[0] = R.dim_[0];
        dim_[1] = R.dim_[1];
        base_   = R.base_;
        nz_     = R.nz_;
        val_    = R.val_;
        rowptr_ = R.rowptr_;
        colind_ = R.colind_;
        return *this;
}

/***************************/
/* newsize()               */
/***************************/

CompRow_Mat_double& CompRow_Mat_double::newsize(int M, int N, int nz)
{
        dim_[0] = M;
        dim_[1] = N;
        nz_     = nz;
        val_.newsize(nz);
        rowptr_.newsize(M+1);
        colind_.newsize(nz);
        return *this;
}

/*********************/
/*   Array access    */
/*********************/


double CompRow_Mat_double::operator()(int i, int j)  const
{
        for (int t=rowptr_(i); t<rowptr_(i+1); t++)
           if (colind_(t) == j) return val_(t);
        if (i < dim_[0] && j < dim_[1]) return 0.0;
        else
        {
            std::cerr << "Array accessing exception -- out of bounds." << "\n";
            return (0);
        }
}


double& CompRow_Mat_double::set(int i, int j)
{
        for (int t=rowptr_(i); t<rowptr_(i+1); t++)
           if (colind_(t) == j) return val_(t);
        std::cerr << "Array element (" << i << "," << j << 
                                ") not in sparse structure -- cannot assign."
             << "\n";          
        exit(1);
    return val_(0);   // // return to suppress compiler warning message
}

/*************/
/*   I/O     */
/*************/
using namespace std;

std::ostream& operator << (std::ostream & os, const CompRow_Mat_double & mat)
{

        int M = mat.dim(0);
        int N = mat.dim(1);
        int rowp1, colp1;
        int flag = 0;


        std::ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        std::ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
        
	int oldp = os.precision(12);


//      Loop through rows...
        for (int i = 0; i < M ; i++)
           for (int j=mat.row_ptr(i);j<mat.row_ptr(i+1);j++)
           {   
              rowp1 =  i + 1;
              colp1 =  mat.col_ind(j) + 1;
              if ( rowp1 == M && colp1 == N ) flag = 1;
              os.width(14);
              os <<  rowp1 ; os << "    " ;
              os.width(14);
              os <<  colp1 ; os << "    " ;
              os.width(20);
              os <<  mat.val(j) << "\n";
           }

        if (flag == 0)
        {
           os.width(14);
           os <<  M ; os << "    " ;
           os.width(14);
           os <<  N ; os << "    " ;
           os.width(20);
           os <<  mat(M-1,N-1) << "\n";
        }

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);

        return os;
}

/***************************************/
/*  Matrix-MV_Vector multiplication    */
/***************************************/

VECTOR_double CompRow_Mat_double::operator*(const VECTOR_double &x) 
        const
{
        int M = dim_[0];
        int N = dim_[1];

//      Check for compatible dimensions:
        if (x.size() != N) 
        {
           std::cerr << "Error in CompCol Matvec -- incompatible dimensions." 
                << "\n";
           exit(1);
           return x;      // return to suppress compiler warning message
         }

         VECTOR_double result(M, 0.0);
         VECTOR_double work(M);
  
         int descra[9];
         descra[0] = 0;
         descra[1] = 0;
         descra[2] = 0;

         F77NAME(dcsrmm) (0, M, 1, N, 1.0,
                  descra, &val_(0), &colind_(0), &rowptr_(0),
                  &x(1), N, 1.0, &result(0), M,
                  &work(0), M);

         return result;
}


/*************************************************/
/* Matrix-Transpose-MV_Vector multiplication...  */
/*************************************************/

VECTOR_double CompRow_Mat_double::trans_mult(const VECTOR_double &x)  
          const
{
          int M = dim_[0];
          int N = dim_[1];

//        Check for compatible dimensions:
          if (x.size() != M) 
          {
             std::cerr << "Error in CompCol Matvec -- incompatible dimensions." 
                  << "\n";
             exit(1);
             return x;   // return to suppress compiler warning message
           }

           VECTOR_double result(N, 0.0);
           VECTOR_double work(N);
  
           int descra[9];
           descra[0] = 0;
           descra[1] = 0;
           descra[2] = 0;

           F77NAME(dcsrmm) (1, N, 1, M, 1.0,
                        descra, &val_(0), &colind_(0), &rowptr_(0),
                    &x(0), M, 1.0, &result(1), N,
                    &work(1), N);

           return result;
}


