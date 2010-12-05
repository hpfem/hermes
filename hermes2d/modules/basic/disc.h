// This is a simple utility that reads strings
// and numbers from a text file. It accepts comments 
// after the '#' symbol.

#include <stdio.h>

// strings:

int Get(FILE *f, char *what);
int Get(FILE *f, char **what);

// first part: 

int Get(FILE *f, short int *what);
int Get(FILE *f, unsigned short int *what);
int Get(FILE *f, int *what);
int Get(FILE *f, unsigned int *what);
int Get(FILE *f, long int *what);
int Get(FILE *f, unsigned long int *what);
int Get(FILE *f, float *what);
int Get(FILE *f, double *what);
int Get(FILE *f, long double *what);
int Get(FILE *f, bool *what);

// second part

int GetFieldLength(FILE *f, int *FieldLength);
void GoToNewOne(FILE *f);

int Get(FILE *f, int *FieldLength, short int **Field);
int Get(FILE *f, int *FieldLength, unsigned short int **Field);
int Get(FILE *f, int *FieldLength, int **Field);
int Get(FILE *f, int *FieldLength, unsigned int **Field);
int Get(FILE *f, int *FieldLength, long int **Field); 
int Get(FILE *f, int *FieldLength, unsigned long int **Field); 
int Get(FILE *f, int *FieldLength, float **Field); 
int Get(FILE *f, int *FieldLength, double **Field); 
int Get(FILE *f, int *FieldLength, long double **Field);

// END OF THE INCLUDE FILE





