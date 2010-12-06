#ifndef _DIRAC_FORMS_H_
#define _DIRAC_FORMS_H_

#include "hermes1d.h"

void assemble_dirac(Mesh *mesh, Matrix *A, Matrix *B, int _kappa, int _Z);

#endif
