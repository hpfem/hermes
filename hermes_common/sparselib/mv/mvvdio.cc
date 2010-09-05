#include <iostream>                                 
#include <stdio.h>
#include "mvvd.h"
#include "mvvi.h"
#include "iotext.h"

int writetxtfile_vec(const char *filename, const MV_Vector_double &A)
{
    FILE *in_file;

    in_file = fopen(filename, "w");
    if (in_file == NULL)
    {
       fprintf(stderr,"Cannot open file:  %s\n", filename );
       exit(1);
    }
    int N = A.size();

    for (int i=0; i<N; i++)
    {
        fprintf(in_file, "%g\n", A(i));
    }

    fclose(in_file);

    return 0;
}

int readtxtfile_vec(const char *filename, MV_Vector_double *Aptr)
{
    FILE *in_file;
    char line[82];
    char *line_ptr;
    int count=0;
    double tmp;
    MV_Vector_double &A = *Aptr;

    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       std::cerr << "Cannot open file: " << filename << "\n";
       exit(1);
    }

    // one vector element per line
    while ( line_ptr = fgets(line, 82, in_file))
        if (sscanf(line_ptr, "%lg", &tmp) >= 1) count++;

    rewind(in_file);

    A.newsize(count);
    for (int i=0; i< count; i++)
    {
        if (fscanf(in_file, "%lg", &A(i)) < 1)
        {
            printf("Error reading %s\n", filename);
            exit(1);
        }
    }

    fclose(in_file);

    return 0;
}


int writetxtfile_vec(const char *filename, const MV_Vector_int &A)
{
    FILE *in_file;

    in_file = fopen(filename, "w");
    if (in_file == NULL)
    {
       fprintf(stderr,"Cannot open file:  %s\n", filename );
       exit(1);
    }
    int N = A.size();

    for (int i=0; i<N; i++)
    {
        fprintf(in_file, "%d\n", A(i));
    }

    fclose(in_file);

    return 0;
}

int readtxtfile_vec(const char *filename, MV_Vector_int *Aptr)
{
    FILE *in_file;
    char line[82];
    char *line_ptr;
    int count=0;
    MV_Vector_int &A = *Aptr;
    int tmp;

    in_file = fopen( filename, "r");
    if (in_file == NULL)
    {
       std::cerr << "Cannot open file: " << filename << "\n";
       exit(1);
    }

    // one vector element per line
    while ( line_ptr = fgets(line, 82, in_file))
        if (sscanf(line_ptr, "%d", &tmp) >= 1) count++;

    rewind(in_file);

    A.newsize(count);
    for (int i=0; i< count; i++)
    {
        if (fscanf(in_file, "%d", &A(i)) < 1)
        {
            printf("Error reading %s\n", filename);
            exit(1);
        }
    }

    fclose(in_file);

    return 0;
}
