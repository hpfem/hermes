#ifndef _SCH_FORMS_H_
#define _SCH_FORMS_H_

#include "hermes1d.h"

double lhs(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

double rhs(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

#define eqn_type_R  0
#define eqn_type_rR 1

void assemble_schroedinger(Mesh *mesh, Matrix *A, Matrix *B, int _l, int Z,
        int equation_type=eqn_type_R);

#endif
