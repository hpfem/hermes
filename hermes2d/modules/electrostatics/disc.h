// last changes 14.6.1995 19:00

#include <stdio.h>

// strings:

int Get(FILE *f, char *what);
int Get(FILE *f, char **what);
void Put(FILE *f, char *what);
void Put(char *what);
void PutNl(FILE *f, char *what);
void PutNl(char *what);

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

// third part:

void Put(FILE *f, short int what);
void Put(FILE *f, unsigned short int what); 
void Put(FILE *f, int what); 
void Put(FILE *f, unsigned int what);
void Put(FILE *f, long int what);
void Put(FILE *f, unsigned long int what);
void Put(FILE *f, float what);
void Put(FILE *f, double what);
void Put(FILE *f, long double what);

// fourth part:

void PutNl(FILE *f, short int what);
void PutNl(FILE *f, unsigned short int what); 
void PutNl(FILE *f, int what);
void PutNl(FILE *f, unsigned int what); 
void PutNl(FILE *f, long int what); 
void PutNl(FILE *f, unsigned long int what);
void PutNl(FILE *f, float what); 
void PutNl(FILE *f, double what);
void PutNl(FILE *f, long double what);

// fifth part:

void Put(short int what); 
void Put(unsigned short int what);
void Put(int what);
void Put(unsigned int what);
void Put(long int what);
void Put(unsigned long int what);
void Put(float what); 
void Put(double what);
void Put(long double what); 

// sixth part:

void PutNl(short int what); 
void PutNl(unsigned short int what); 
void PutNl(int what);
void PutNl(unsigned int what); 
void PutNl(long int what);
void PutNl(unsigned long int what);
void PutNl(float what);
void PutNl(double what);
void PutNl(long double what);

// END OF THE INCLUDE FILE





