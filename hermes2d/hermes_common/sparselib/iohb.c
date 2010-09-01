/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Update: (10/28/2008): removed <malloc.h> dependency.                    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*      I/O for Harwell-Boeing files                                        */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                          */
/*  Several user I/O functions for Harwell-Boeing files are provided,       */
/*  together with various utility functions.                                */
/*  For a description of the Harwell-Boeing standard, see:                  */
/*                                                                          */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Function:                                                               */
/*  void readHB_newmat_double(char *filename, int *M, int *N, *int nz,      */
/*        int **colptr, int **rowind,  double**val)                         */
/*                                                                          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Function:                                                               */
/*  void readHB_info(char *filename, int *M, int *N, int *nz, int *nrhs)    */
/*                                                                          */
/*  Description:                                                            */
/*                                                                          */
/*  The readHB_info function opens and reads the header information from    */
/*  the specified Harwell-Boeing file, and reports back the number of rows  */
/*  and columns in the stored matrix (M and N), the number of nonzeros in   */
/*  the matrix (nz), and the number of right-hand-sides stored along with   */
/*  the matrix (nrhs).                                                      */
/*  The optional verbose parameter, if set to 1,  triggers the output of    */
/*  more detailed information from the header, including title, etc.        */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Function:                                                               */
/*  void readHB_mat_double(char *filename, int *colptr, int *rowind,        */
/*             double*val)                                                  */
/*                                                                          */
/*  Description:                                                            */
/*                                                                          */
/*  This function opens and reads the specified file, interpreting its      */
/*  contents as a sparse matrix stored in the Harwell/Boeing standard       */
/*  format and creating compressed column storage scheme vectors to hold    */
/*  the index and nonzero value information.                                */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Function:                                                               */
/*  void readHB_rhs_double(char *filename, double *b, int j)                */
/*                                                                          */
/*  Description:                                                            */
/*                                                                          */
/*  This function opens and reads the specified file, returning a right-    */
/*  hand-side vector b.  If the file specifies multiple right-hand-sides    */
/*  (that is, a right-hand-side matrix B), then the optional argument j     */
/*  may be used to select the (j+1)st column of B.                          */
/*                                                                          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Function:                                                               */
/*  void writeHB_mat_double(char *filename, int M, int N, int nz,           */
/*         *colptr,  int *rowind, double *val, int nrhs, double *rhs,       */
/*         char *Title, char *Key)                                          */
/*                                                                          */
/*  Description:                                                            */
/*                                                                          */
/*  The writeHB function opens the named file and writes the specified      */
/*  matrix and optional right-hand-side(s) to that file in Harwell-Boeing   */
/*  format.                                                                 */
/*                                                                          */
/****************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void nullchk();
void ParseIfmt();
void ParseRfmt();
void convertDtoE();
void readHB_header();
void readHB_mat_double();
void readHB_rhs_double();
void writeHB_mat_double();
char * substr();
char * substr_after();
char * substr_before();
char * substr_through();
void upcase();
int my_index();

void readHB_info(filename, M, N, nz, nrhs)
char *filename; int *M; int *N; int *nz; int *nrhs; 
{
/****************************************************************************/
/*  The readHB_info function opens and reads the header information from    */
/*  the specified Harwell-Boeing file, and reports back the number of rows  */
/*  and columns in the stored matrix (M and N), the number of nonzeros in   */
/*  the matrix (nz), and the number of right-hand-sides stored along with   */
/*  the matrix (nrhs).                                                      */
/*  The optional verbose parameter, if set to 1,  triggers the output of    */
/*  more detailed information from the header, including title, etc.        */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    int verbose = 0;
    FILE *in_file;
    char Title[73], Key[9], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd; 
    int Nrow, Ncol, Nnzero;
    int Nrhs;
    
    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       printf("Cannot open file: %s\n",filename);
       exit(1);
    }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt, 
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

    *M    = Nrow;
    *N    = Ncol;
    *nz   = Nnzero;
    if (Rhscrd != 0) {*nrhs = Nrhs;}
                else {*nrhs = 0;}

/*  In verbose mode, print some of the header information:   */
    if (verbose == 1)
    {
        printf("Reading from Harwell-Boeing file %s...\n",filename);
        printf("Title: %s \n",Title);
        printf("Key:   %s \n",Key);
        printf("The stored matrix is %i by %i with %i nonzeros.\n", 
                *M, *N, *nz );
        printf("%i right-hand--side(s) are stored.\n",*nrhs);
    }
 
    return;

}


/*************************************************************************/
/*  Read header information from the named H/B file...                   */
/*************************************************************************/

void  readHB_header(in_file, Title, Key, Type, 
                    Nrow, Ncol, Nnzero, Nrhs,
                    Ptrfmt, Indfmt, Valfmt, Rhsfmt, 
                    Ptrcrd, Indcrd, Valcrd, Rhscrd, 
                    Rhstype)
FILE *in_file; char *Title; char *Key; char *Type; 
int *Nrow; int *Ncol; int *Nnzero; int *Nrhs;
char *Ptrfmt; char *Indfmt; char *Valfmt; char *Rhsfmt;
int *Ptrcrd; int *Indcrd; int *Valcrd; int *Rhscrd;
char *Rhstype;
{
    char line[82];
    char *line_ptr;
    int Totcrd;
    int Neltvl, Nrhsix;

/*  First line:   */
    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    (void) sscanf(line, "%72c %8c", Title, Key);

/*  Second line:  */
    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    (void) sscanf(line, "%i %i %i %i %i", &Totcrd, Ptrcrd, Indcrd, 
                                          Valcrd, Rhscrd);

/*  Third line:   */
    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    (void) sscanf(line, "%3c %i %i %i %i", Type, Nrow, Ncol, 
                                          Nnzero, &Neltvl);

/*  Fourth line:  */
    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    (void) sscanf(line, "%s %s %s %s", Ptrfmt, Indfmt, Valfmt, Rhsfmt);
   
/*  (Optional) Fifth line: */
    if (*Rhscrd != 0 )
    { 
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       (void) sscanf(line, "%3c %i %i", Rhstype, Nrhs, &Nrhsix);
    }

    return;
}



void readHB_mat_double(filename, colptr, rowind, val)
char *filename; int colptr[]; int rowind[]; double val[];
{
/****************************************************************************/
/*  This function opens and reads the specified file, interpreting its      */
/*  contents as a sparse matrix stored in the Harwell/Boeing standard       */
/*  format and creating compressed column storage scheme vectors to hold    */
/*  the index and nonzero value information.                                */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    FILE *in_file;
    char Title[73], Key[9], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero;
    int Nrhs;
    char line[82];
    char* line_ptr;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline, Valwidth;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    int i, ind, col;
    int count;
    char* ThisElement;

    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       printf("Cannot open file: %s\n",filename);
       exit(1);
    }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

/*  Parse the array input formats from Line 3 of HB file  */
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valflag);

/*  Read column pointer array:   */

    count=0;
    for (i=0;i<Ptrcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       col =  0;
       for (ind = 0;ind<Ptrperline;ind++)
       {
          if (count > Ncol) break;
          ThisElement = substr(line,col,Ptrwidth);
          colptr[count] = atoi(ThisElement)-1;
          count++; col += Ptrwidth;
       }
    }

/*  Read row index array:  */

    count = 0;
    for (i=0;i<Indcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       col =  0;
       for (ind = 0;ind<Indperline;ind++)
       {
          if (count == Nnzero) break;
          ThisElement = substr(line,col,Indwidth);
          rowind[count] = atoi(ThisElement)-1;
          count++; col += Indwidth;
       }
    }

/*  Read array of values:  */

    count = 0;
    for (i=0;i<Valcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       if (Valflag == 'D') convertDtoE(line);
       col =  0;
       for (ind = 0;ind<Valperline;ind++)
       {
          if (count == Nnzero) break;
          ThisElement = substr(line,col,Valwidth);
          val[count] = atof(ThisElement);
          count++; col += Valwidth;
       }
    }

    return;
}


void readHB_mat_float(filename, colptr, rowind, val)
char *filename; int *colptr; int *rowind; float *val;
{
/****************************************************************************/
/*  This function opens and reads the specified file, interpreting its      */
/*  contents as a sparse matrix stored in the Harwell/Boeing standard       */
/*  format and creating compressed column storage scheme vectors to hold    */
/*  the index and nonzero value information.                                */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    FILE *in_file;
    char Title[73], Key[9], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero;
    int Nrhs;
    char line[82];
    char* line_ptr;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline, Valwidth;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    int i, ind, col;
    int count;
    char* ThisElement;

    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       printf("Cannot open file: %s\n",filename);
       exit(1);
    }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

/*  Parse the array input formats from Line 3 of HB file  */
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valflag);

/*  Read column pointer array:   */

    count=0;
    for (i=0;i<Ptrcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       col =  0;
       for (ind = 0;ind<Ptrperline;ind++)
       {
          if (count > Ncol) break;
          ThisElement = substr(line,col,Ptrwidth);
          colptr[count] = atoi(ThisElement)-1;
          count++; col += Ptrwidth;
       }
    }

/*  Read row index array:  */

    count = 0;
    for (i=0;i<Indcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       col =  0;
       for (ind = 0;ind<Indperline;ind++)
       {
          if (count == Nnzero) break;
          ThisElement = substr(line,col,Indwidth);
          rowind[count] = atoi(ThisElement)-1;
          count++; col += Indwidth;
       }
    }

/*  Read array of values:  */

    count = 0;
    for (i=0;i<Valcrd;i++)
    {
       line_ptr = fgets(line, 82, in_file);
       nullchk(line_ptr);
       if (Valflag == 'D') convertDtoE(line);
       col =  0;
       for (ind = 0;ind<Valperline;ind++)
       {
          if (count == Nnzero) break;
          ThisElement = substr(line,col,Valwidth);
          val[count] = atof(ThisElement);
          count++; col += Valwidth;
       }
    }

    return;
}


void readHB_rhs_double(filename, b, j)
char *filename; double b[]; int j;
{
/****************************************************************************/
/*  This function opens and reads the specified file, returning a right-    */
/*  hand-side vector b.  If the file specifies multiple right-hand-sides    */
/*  (that is, a right-hand-side matrix B), then the optional argument j     */
/*  may be used to select the (j+1)st column of B.                          */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    int i, n, null_entries, lines_left;
    FILE *in_file;
    char Title[73], Key[9], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero;
    int Nrhs;
    int Rhsperline, Rhswidth;
    int Rhsflag;
    char buffer[BUFSIZ];
    int ind, col;
    int count;
    char *line_ptr;
    char line[80];
    char *ThisElement;

    if ((in_file = fopen( filename, "r")) == NULL) 
    {
      printf("Cannot open file: %s\n",filename);
      exit(1);
     }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

    if (Nrhs <= 0)
    {
      printf("Attempt to read rhs when none is present.\n");
      exit(1);
      return;
    }
    if (Rhstype[0] != 'F')
    {
      printf("Attempt to read rhs which is not stored in Full form.\n");
      printf("Rhs must be specified as full. \n");
      exit(1);
      return;
    }

    ParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsflag);

/*  Lines to skip before starting to read RHS values... */
    n = Ptrcrd + Indcrd + Valcrd + j*Nrow/Rhsperline;

/*  number of entries on the line of interest to skip before */
/*  reading RHS values...                                    */
    null_entries = j*Nrow%Rhsperline;
    lines_left = (int) ( .5 + (double) (Nrow-Rhsperline+null_entries)/ 
                                            (double) Rhsperline );

    for (i = 0; i < n; i++)
      fgets(buffer, BUFSIZ, in_file);

    count = 0;

/*  Handle first line separately, in case j != 0 and the right-hand-side */
/*  of interest starts mid-line (after null_entries entries)...          */

    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    col = 0;
    for (ind = 0; ind < Rhsperline; ind++)
    {
      if (ind < null_entries)
        col += Rhswidth;
      else
      {
        if (count > Nrow-1)
          break;
        ThisElement = substr(line, col, Rhswidth);
        b[count] = atof(ThisElement);
        count++; col += Rhswidth;
      }
    }

/*  Handle subsequent lines...              */
    for (i = 0; i < lines_left; i++)
    {
      line_ptr = fgets(line, 82, in_file);
      nullchk(line_ptr);
      col = 0;
      for (ind = 0; ind < Rhsperline; ind++)
      {
         if (count > Nrow-1)
           break;
         ThisElement = substr(line, col, Rhswidth);
         b[count] = atof(ThisElement);
         count++; col += Rhswidth;
      }
    }

    return;
}

void readHB_rhs_float(filename, b, j)
char *filename; float b[]; int j;
{
/****************************************************************************/
/*  This function opens and reads the specified file, returning a right-    */
/*  hand-side vector b.  If the file specifies multiple right-hand-sides    */
/*  (that is, a right-hand-side matrix B), then the optional argument j     */
/*  may be used to select the (j+1)st column of B.                          */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    int i, n, null_entries, lines_left;
    FILE *in_file;
    char Title[73], Key[9], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero;
    int Nrhs;
    int Rhsperline, Rhswidth;
    int Rhsflag;
    char buffer[BUFSIZ];
    int ind, col;
    int count;
    char *line_ptr;
    char line[80];
    char *ThisElement;

    if ((in_file = fopen( filename, "r")) == NULL) 
    {
      printf("Cannot open file: %s\n",filename);
      exit(1);
     }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

    if (Nrhs <= 0)
    {
      printf("Attempt to read rhs when none is present.\n");
      exit(1);
      return;
    }
    if (Rhstype[0] != 'F')
    {
      printf("Attempt to read rhs which is not stored in Full form.\n");
      printf("Rhs must be specified as full. \n");
      exit(1);
      return;
    }

    ParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsflag);

/*  Lines to skip before starting to read RHS values... */
    n = Ptrcrd + Indcrd + Valcrd + j*Nrow/Rhsperline;

/*  number of entries on the line of interest to skip before */
/*  reading RHS values...                                    */
    null_entries = j*Nrow%Rhsperline;
    lines_left = (int)( .5 + (double) (Nrow-Rhsperline+null_entries)/
                                             (double) Rhsperline);

    for (i = 0; i < n; i++)
      fgets(buffer, BUFSIZ, in_file);

    count = 0;

/*  Handle first line separately, in case j != 0 and the right-hand-side */
/*  of interest starts mid-line (after null_entries entries)...          */

    line_ptr = fgets(line, 82, in_file);
    nullchk(line_ptr);
    col = 0;
    for (ind = 0; ind < Rhsperline; ind++)
    {
      if (ind < null_entries)
        col += Rhswidth;
      else
      {
        if (count > Nrow-1)
          break;
        ThisElement = substr(line, col, Rhswidth);
        b[count] = atof(ThisElement);
        count++; col += Rhswidth;
      }
    }

/*  Handle subsequent lines...              */
    for (i = 0; i < lines_left; i++)
    {
      line_ptr = fgets(line, 82, in_file);
      nullchk(line_ptr);
      col = 0;
      for (ind = 0; ind < Rhsperline; ind++)
      {
         if (count > Nrow-1)
           break;
         ThisElement = substr(line, col, Rhswidth);
         b[count] = atof(ThisElement);
         count++; col += Rhswidth;
      }
    }

    return;
}


void writeHB_mat_double(filename, M, N, nz, colptr, rowind, val, nrhs, rhs, 
                                                                  Title, Key)
char *filename; int M, N, nz, colptr[], rowind[]; double val[]; int nrhs;
double rhs[]; char *Title; char *Key;
{
/****************************************************************************/
/*  The writeHB function opens the named file and writes the specified      */
/*  matrix and optional right-hand-side(s) to that file in Harwell-Boeing   */
/*  format.                                                                 */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/****************************************************************************/


    FILE *out_file;

    int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    char *ptrfmt, *indfmt, *valfmt, *rhsfmt;
    char *Type;
    char *Titlefill, *Keyfill;
    int filllen, i;
    int entry, finfo;

    out_file = fopen( filename, "w");

    ptrcrd = (N+1)/8;
    if ( (N+1)%8 != 0) ptrcrd++;

    indcrd = nz/8;
    if ( nz%8 != 0) indcrd++;

    valcrd = nz/4;
    if ( nz%4 != 0) valcrd++;

    rhscrd = nrhs*M/4; 
    if ( nrhs*M%4 != 0) rhscrd++;

    totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;

    ptrfmt = "(8I10)          ";
    indfmt = ptrfmt;
    valfmt = "(4E20.16)           ";
    rhsfmt = valfmt;

    Type = "RUA";

/*  Print header information:  */

    fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
            ptrcrd, indcrd, valcrd, rhscrd);
    fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
    fprintf(out_file,"%16s%16s%20s%20s\n", ptrfmt, indfmt, valfmt, rhsfmt);
    if ( nrhs != 0 ) {
/*     Print optional fifth header line for right-hand-side information : */
       fprintf(out_file,"F             %d\n", nrhs);
    }

/*  Print column pointers:   */
    for (i=0;i<N+1;i++)
    {
       entry = colptr[i]+1;
       fprintf(out_file,"%10d",entry);
       if ( (i+1)%8 == 0 ) fprintf(out_file,"\n");
    }

   if ( (N+1) % 8 != 0 ) fprintf(out_file,"\n");

/*  Print row indices:       */
    for (i=0;i<nz;i++)
    {
       entry = rowind[i]+1;
       fprintf(out_file,"%10d",entry);
       if ( (i+1)%8 == 0 ) fprintf(out_file,"\n");
    }

   if ( nz % 8 != 0 ) fprintf(out_file,"\n");

/*  Print values:            */
    for (i=0;i<nz;i++)
    {
       fprintf(out_file,"% 20.12E",val[i]);
       if ( (i+1)%4 == 0 ) fprintf(out_file,"\n");
    }

    if ( nz % 4 != 0 ) fprintf(out_file,"\n");

/*  Print right hand sides:  */
    if ( nrhs > 0 ) {
       for (i=0;i<nrhs*M;i++)
       {
          fprintf(out_file,"% 20.12E",rhs[i]);
          if ( (i+1)%4 == 0 ) fprintf(out_file,"\n");
       }
    }

    finfo = fclose(out_file);
    if (finfo != 0) printf("Error closing file in writeHB_mat_double().\n");

    return;
}


void writeHB_mat_float(filename, M, N, nz, colptr, rowind, val, nrhs, rhs, 
                                                                  Title, Key)
char *filename; int M, N, nz, colptr[], rowind[]; float val[]; int nrhs;
float rhs[]; char *Title; char *Key;
{
/****************************************************************************/
/*  The writeHB function opens the named file and writes the specified      */
/*  matrix and optional right-hand-side(s) to that file in Harwell-Boeing   */
/*  format.                                                                 */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/****************************************************************************/


    FILE *out_file;

    int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    char *ptrfmt, *indfmt, *valfmt;
    char *Type;
    char *Titlefill, *Keyfill;
    int filllen, i;
    int entry, finfo;

    out_file = fopen( filename, "w");

    ptrcrd = (N+1)/8;
    if ( (N+1)%8 != 0) ptrcrd++;

    indcrd = nz/8;
    if ( nz%8 != 0) indcrd++;

    valcrd = nz/4;
    if ( nz%4 != 0) valcrd++;

    rhscrd = nrhs*M/4; 
    if ( nrhs*M%4 != 0) rhscrd++;

    totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;

    ptrfmt = "(8I10)          ";
    indfmt = ptrfmt;
    valfmt = "(4E20.16)       ";

    Type = "RUA";

/*  Print header information:  */

    fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
            ptrcrd, indcrd, valcrd, rhscrd);
    fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
    fprintf(out_file,"%16s%16s%16s\n", ptrfmt, indfmt, valfmt);

/*  Print column pointers:   */
    for (i=0;i<N+1;i++)
    {
       entry = colptr[i]+1;
       fprintf(out_file,"%10d",entry);
       if ( (i+1)%8 == 0 ) fprintf(out_file,"\n");
    }

   if ( (N+1) % 8 != 0 ) fprintf(out_file,"\n");

/*  Print row indices:       */
    for (i=0;i<nz;i++)
    {
       entry = rowind[i]+1;
       fprintf(out_file,"%10d",entry);
       if ( (i+1)%8 == 0 ) fprintf(out_file,"\n");
    }

   if ( nz % 8 != 0 ) fprintf(out_file,"\n");

/*  Print values:            */
    for (i=0;i<nz;i++)
    {
       fprintf(out_file,"% 20.12E",val[i]);
       if ( (i+1)%4 == 0 ) fprintf(out_file,"\n");
    }

    if ( nz % 4 != 0 ) fprintf(out_file,"\n");

/*  Print right hand sides:  */
    if ( nrhs > 0 ) {
       for (i=0;i<nrhs*M;i++)
       {
          fprintf(out_file,"% 20.12E",rhs[i]);
          if ( (i+1)%4 == 0 ) fprintf(out_file,"\n");
       }
    }

    finfo = fclose(out_file);
    if (finfo != 0) printf("Error closing file in writeHB_mat_float().\n");

    return;
}

/****************************************************************************/

              /********************************************/
              /*  Utility functions supporting readHB:    */
              /********************************************/

/****************************************************************************/

/**************************************************/
/*  Check to ensure header information is present */
/**************************************************/
void nullchk(line_ptr)
char *line_ptr;
{
    if (line_ptr == NULL)
    {
          printf("Cannot complete reading file information.\n ");
          exit(1);
    }
}

/*************************************************/
/*  Parse an *integer* format field to determine */
/*  width and number of elements per line.       */
/*************************************************/
void ParseIfmt(fmt, perline, width)
char *fmt; int *perline; int *width;
{
    char *PPL, *WID;
    upcase(fmt);
    PPL = substr_before(fmt,'I');
    PPL = substr(PPL,1,strlen(PPL)-1);
    *perline = atoi(PPL);
    WID = substr_after(fmt,'I');
    WID = substr(WID,0,strlen(WID)-1);
    *width = atoi(WID);
}

/*************************************************/
/*  Parse an *real* format field to determine    */
/*  width and number of elements per line.       */
/*  Also sets flag indicating 'E' or 'D' format. */
/*************************************************/
void ParseRfmt(fmt, perline, width, flag)
char *fmt; int *perline; int *width; int *flag;
{
    int foundE, foundD, foundF, foundP;
    char *PPL, *WID;

    if (fmt == NULL )
    {
      *perline = 0; *width = 0; *flag = '\0'; return;
    }

    upcase(fmt);
    foundP = my_index(fmt,'P');
    foundE = my_index(fmt,'E');
    foundD = my_index(fmt,'D');
    foundF = my_index(fmt,'F');      /* Fixed format */
    if (foundP != -1 )            /* Remove any scaling factor, which */
    {                             /* affects output only, not input */
      fmt = fmt + foundP + 1;
    }
    if (foundE != -1 )
    {
      *flag = 'E';
      PPL = substr_before(fmt,'E');
      PPL = substr(PPL,1,strlen(PPL)-1);
      *perline = atoi(PPL);
      WID = substr_after(fmt,'E');
      WID = substr_through(WID,'.');
      WID = substr(WID,0,strlen(WID)-1);
      *width = atoi(WID);
    }
    else if (foundD != -1)
    {
      *flag = 'D';
      PPL = substr_before(fmt,'D');
      PPL = substr(PPL,1,strlen(PPL)-1);
      *perline = atoi(PPL);
      WID = substr_after(fmt,'D');
      WID = substr_through(WID,'.');
      WID = substr(WID,0,strlen(WID)-1);
      *width = atoi(WID);
    }
    else if (foundF != -1)
    {
      *flag = 'F';
      PPL = substr_before(fmt,'F');
      PPL = substr(PPL,1,strlen(PPL)-1);
      *perline = atoi(PPL);
      WID = substr_after(fmt,'F');
      WID = substr_through(WID,'.');
      WID = substr(WID,0,strlen(WID)-1);
      *width = atoi(WID);
    }
    else
    {
      printf("Real format in H/B file not supported.\n");
      exit(1);
    }

}


/*****************************************************/
/*  Converts a line with real data in 'D' format     */
/*  to 'E' format by simply changing the character   */
/*  'D' to 'E' by adding 1 to it.                    */
/*****************************************************/
void convertDtoE(line)
char *line;
{
    int len, i;
    len = strlen(line);
    for (i=0;i<len;i++)
       if ( line[i] == 'D' || line[i] == 'd' ) line[i] = line[i]+1;
    return;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*       small string processing library for readhb()            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*****************************************************************/
/*  This file contains some minimal string functions to support  */
/*  support parsing the data formats in a Harwell-Boeing file.   */
/*****************************************************************/

char* substr(S, pos, len)
char *S;
int pos, len;
{
    int i;
    char *SubS;
    SubS = malloc(len+1);
    if ( (pos+len) > strlen(S))
    {
       printf("Error: Substring (%s, %d, %d) will read beyond string boundary.\n",
               S, pos, len);
       exit(1);
    }
    for (i=0;i<len;i++)
       SubS[i] = S[pos+i];
    SubS[len] = '\0';
    return SubS;
}

char* substr_after(S, M)
char *S; char M;
{
/*  Return the substring of S from the char after the first   */
/*  occurrence of M to the end of the string                  */
    int i, pos, SubSlen;
    char * SubS;
    pos = 0;
    for (i=0;i< strlen(S);i++)
    {
       if ( S[i] == M )
       {
          pos = i+1;
          break;
       }
    }
    if (pos != 0)
    {
       SubSlen = strlen(S) - pos;
       SubS = malloc(SubSlen);
       for (i=0;i< SubSlen;i++)
          SubS[i] = S[pos+i];
       SubS[SubSlen] = '\0';
       return SubS;
     }
     else
     {
       printf("Character %i not found in input string.\n", M);
       exit(1);
     }

       return NULL;       
 }



char* substr_before(S, M)
char *S; char M;
{
/*  Return the substring of S from the first char of S to the   */
/*  char before the first occurrence of M                       */
    int i, pos, SubSlen;
    char* SubS;
    pos = 0;
    for (i=0;i< strlen(S);i++)
    {
       if ( S[i] == M )
       {
          pos = i-1;
          break;
       }
    }
    if (pos != 0)
    {
       SubSlen = pos+1;
       SubS = malloc(SubSlen);
       for (i=0;i< SubSlen;i++)
          SubS[i] = S[i];
       SubS[SubSlen] = '\0';
       return SubS;
     }
     else
     {
       printf("Character %i not found in input string.\n", M);
       exit(1);
     }
     return NULL;   
 }


char* substr_through(S, M)
char *S; char M;
{
/*  Similer to substr_before, but include M         */
    int i, pos, SubSlen;
    char *SubS;
    pos = 0;
    for (i=0;i< strlen(S);i++)
    {
       if ( S[i] == M )
       {
          pos = i;
          break;
       }
    }
    if (pos != 0)
    {
       SubSlen = pos+1;
       SubS = malloc(SubSlen);
       for (i=0;i< SubSlen;i++)
          SubS[i] = S[i];
       SubS[SubSlen] = '\0';
       return SubS;
     }
     else
     {
       printf("Character %i not found in input string.\n", M);
       exit(1);
     }
       return NULL;
 }


void upcase(S)
char *S;
{
/*  Convert S to uppercase     */
    int i;
    for (i=0;i< strlen(S);i++)
       if ( S[i] >= 'a' && S[i] <= 'z' ) S[i] = S[i] - 32;
}


int my_index(S, M)
char *S; char M;
{
/*  Return the index of the first occurrence of M in S  */
/*  Return -1 if M does not occur in S                  */
    int i, pos;
    pos = -1;
    for (i=0;i< strlen(S);i++)
    {
       if ( S[i] == M )
       {
          pos = i;
          break;
       }
    }
    return pos;
}


void readHB_newmat_double(filename, M, N, nonzeros, colptr, rowind, val)
char *filename; int *M; int *N; int *nonzeros; int **colptr; 
int **rowind; double **val;
{
    int nrhs;

    readHB_info(filename, M, N, nonzeros, &nrhs);
   *colptr = (int *)malloc((*N+1)*sizeof(int));
   *rowind = (int *)malloc(*nonzeros*sizeof(int));
   *val = (double *)malloc(*nonzeros*sizeof(double));

    readHB_mat_double(filename, *colptr, *rowind, *val);

}

void readHB_newmat_float(filename, M, N, nonzeros, colptr, rowind, val)
char *filename; int *M; int *N; int *nonzeros; int **colptr; 
int **rowind; float **val;
{
    int nrhs;

    readHB_info(filename, M, N, nonzeros, &nrhs);
    *colptr = (int *)malloc((*N+1)*sizeof(int));
    *rowind = (int *)malloc(*nonzeros*sizeof(int));
    *val = (float *)malloc(*nonzeros*sizeof(float));

    readHB_mat_float(filename, *colptr, *rowind, *val);

}




