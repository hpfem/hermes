// last changes 14.6.1995 10:00

# include <stdio.h>
# include <stdlib.h>

// strings:

int Get(FILE *f, char *what) {
  unsigned int i=0;
  int t;
  t = fgetc(f);
  if(t == EOF) return 0;
  while((t == ' ') || (t == '\n') || (t == '\t') ||
        (t == '\r') || (t == '#')) {
    while((t == ' ') || (t == '\n') || (t == '\t') || (t == '\r')) {
      t = fgetc(f);
      if(t == EOF) return 0;
    };
    if(t == '#') {
      do {
        t = fgetc(f);
        if(t == EOF) return 0;
      } while(t != '\n');
    };
  };
  do {
    what[i++] = t;
    t = fgetc(f);
  }
  while (!((t == ' ') || (t == '\n') || (t == '\t') ||
          (t == '\r') || (t == EOF)));
  if(t == '\n') fseek(f, -1, SEEK_CUR);
  what[i] = '\0';
  return 1;
}

int Get(FILE *f, char **what) {
  *what = (char*)malloc(255);
  if(!Get(f, *what)) return 0;
  return 1;
}

void Put(FILE *f, char *what) {
  fprintf(f, "%s ", what);
}

void Put(char *what) {
  fprintf(stderr, "%s ", what);
}

void PutNl(FILE *f, char *what) {
  fprintf(f, "%s\n", what);
}

void PutNl(char *what) {
  fprintf(stderr, "%s\n", what);
}

// numbers:

// first part: 
// int Get(FILE *f, type *what);

int Get(FILE *f, short int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (short int)atoi(str);
  return 1;
}

int Get(FILE *f, unsigned short int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (unsigned short int)atoi(str);
  return 1;
}

int Get(FILE *f, int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = atoi(str);
  return 1;
}

int Get(FILE *f, unsigned int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (unsigned int)atoi(str);
  return 1;
}

int Get(FILE *f, long int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = atol(str);
  return 1;
}

int Get(FILE *f, unsigned long int *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (unsigned long int)atol(str);
  return 1;
}

int Get(FILE *f, float *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = atof(str);
  return 1;
}

int Get(FILE *f, double *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (double)atof(str);
  return 1;
}

int Get(FILE *f, long double *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  *what = (long double)atof(str);
  return 1;
}

// second part
// int Get(FILE *f, int FieldLength, type **Field);

void GoToNewOne(FILE *f) {
  int t;
  t = fgetc(f);
  while((t == ' ') || (t == '\t') || (t == '\r') || (t == '#')) {
    if(t == '#') {
      do {
        t = fgetc(f);
        if(t == EOF) return;
      } while(t != '\n');
      break;
    }
    t = fgetc(f);
  }
  fseek(f, -1, SEEK_CUR);
}

int GetFieldLength(FILE *f, int *FieldLength) {
  int t; char str[255];
  unsigned long int Position;
  *FieldLength = 0;
  Position = ftell(f);
  do {
    if(!Get(f, str)) return 0;
    GoToNewOne(f);
    (*FieldLength)++;
    t = fgetc(f);
    if(t != EOF) fseek(f, -1, SEEK_CUR);
  } while(t != EOF && t != '\n');
  fseek(f, Position, SEEK_SET);
  return 1;
}

int Get(FILE *f, int *FieldLength, short int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (short int*) malloc(*FieldLength*sizeof(short int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, unsigned short int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (unsigned short int*) malloc(*FieldLength*sizeof(unsigned short int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (int*) malloc(*FieldLength*sizeof(int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, unsigned int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (unsigned int*) malloc(*FieldLength*sizeof(unsigned int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, long int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (long int*) malloc(*FieldLength*sizeof(long int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, unsigned long int **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (unsigned long int*) malloc(*FieldLength*sizeof(unsigned long int));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, float **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (float*) malloc(*FieldLength*sizeof(float));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, double **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (double*) malloc(*FieldLength*sizeof(double));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

int Get(FILE *f, int *FieldLength, long double **Field) {
  if(!GetFieldLength(f, FieldLength)) return 0;
  *Field = (long double*) malloc(*FieldLength*sizeof(long double));
  if(*Field == NULL) return 0;
  for(int i=0; i<*FieldLength; i++) Get(f, *Field + i);
  return 1;
}

// third part:
// void Put(FILE *f, type what);

void Put(FILE *f, short int what) {
  fprintf(f, "%d ", (int)what);
}

void Put(FILE *f, unsigned short int what) {
  fprintf(f, "%u ", (unsigned int)what);
}

void Put(FILE *f, int what) {
  fprintf(f, "%d ", what);
}

void Put(FILE *f, unsigned int what) {
  fprintf(f, "%u ", what);
}

void Put(FILE *f, long int what) {
  fprintf(f, "%ld ", what);
}

void Put(FILE *f, unsigned long int what) {
  fprintf(f, "%lu ", what);
}

void Put(FILE *f, float what) {
  fprintf(f, "%g ", what);
}

void Put(FILE *f, double what) {
  fprintf(f, "%g ", what);
}

// fourth part:
// void PutNl(FILE *f, type what);

void PutNl(FILE *f, short int what) {
  fprintf(f, "%d\n", (int)what);
}

void PutNl(FILE *f, unsigned short int what) {
  fprintf(f, "%u\n", (unsigned int)what);
}

void PutNl(FILE *f, int what) {
  fprintf(f, "%d\n", what);
}

void PutNl(FILE *f, unsigned int what) {
  fprintf(f, "%u\n", what);
}

void PutNl(FILE *f, long int what) {
  fprintf(f, "%ld\n", what);
}

void PutNl(FILE *f, unsigned long int what) {
  fprintf(f, "%lu\n", what);
}

void PutNl(FILE *f, float what) {
  fprintf(f, "%g\n", what);
}

void PutNl(FILE *f, double what) {
  fprintf(f, "%g\n", what);
}

// fifth part:
// void Put(type what);


void Put(unsigned short int what) {
  fprintf(stderr, "%u ", (unsigned int)what);
}

void Put(int what) {
  fprintf(stderr, "%d ", what);
}

void Put(unsigned int what) {
  fprintf(stderr, "%u ", what);
}

void Put(long int what) {
  fprintf(stderr, "%ld ", what);
}

void Put(unsigned long int what) {
  fprintf(stderr, "%lu ", what);
}

void Put(float what) {
  fprintf(stderr, "%g ", what);
}

void Put(double what) {
  fprintf(stderr, "%g ", what);
}

// sixth part:
// void PutNl(type what);

void PutNl(short int what) {
  fprintf(stderr, "%d\n", (int)what);
}

void PutNl(unsigned short int what) {
  fprintf(stderr, "%u\n", (unsigned int)what);
}

void PutNl(int what) {
  fprintf(stderr, "%d\n", what);
}

void PutNl(unsigned int what) {
  fprintf(stderr, "%u\n", what);
}

void PutNl(long int what) {
  fprintf(stderr, "%ld\n", what);
}

void PutNl(unsigned long int what) {
  fprintf(stderr, "%lu\n", what);
}

void PutNl(float what) {
  fprintf(stderr, "%g\n", what);
}

void PutNl(double what) {
  fprintf(stderr, "%g\n", what);
}

// END OF THE INCLUDE FILE
