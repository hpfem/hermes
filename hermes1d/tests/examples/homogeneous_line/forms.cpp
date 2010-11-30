double jacobian_1_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_1_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_2_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*L*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_2_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -R*u[i]*v[i]*weights[i];
    }
    return val;
};


double jacobian_3_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += +omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_3_3(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};



double jacobian_4_1(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -omega*C*u[i]*v[i]*weights[i];
    }
    return val;
};

double jacobian_4_2(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += -G*weights[i];
    }
    return val;
};

double jacobian_4_4(int num, double *x, double *weights,
                    double *u, double *dudx, double *v, double *dvdx,
                    double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                    void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        val += dudx[i]*v[i]*weights[i];
    }
    return val;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double residual_1(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    int si = 0;      // solution index (only 0 is relevant for this example)
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[si][0][i] - R*u_prev[si][2][i] +L*omega*u_prev[si][3][i])*v[i]*weights[i];
    }
    return val;
};

double residual_2(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    int si = 0;      // solution index (only 0 is relevant for this example)
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[si][1][i] - R*u_prev[si][3][i]-omega*L*u_prev[si][2][i])*v[i]*weights[i];
    }
    return val;
};

double residual_3(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    int si = 0;      // solution index (only 0 is relevant for this example)
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[si][2][i] - G*u_prev[si][0][i] + omega*C*u_prev[si][1][i])*v[i]*weights[i];
    }
    return val;
};


double residual_4(int num, double *x, double *weights,
                  double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                  double *v, double *dvdx, void *user_data)
{
    double val = 0;
    int si = 0;      // solution index (only 0 is relevant for this example)
    for(int i = 0; i<num; i++) {
        val += (du_prevdx[si][3][i] - G*u_prev[si][1][i]-omega*C*u_prev[si][0][i])*v[i]*weights[i];
    }
    return val;
};


double jacobian_surf_right_U_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{     
    return 1*u*v;
}


double jacobian_surf_right_U_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}

double jacobian_surf_right_I_Re(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
    return 1*u*v;
}


double jacobian_surf_right_I_Im(double x, double u, double dudx,
                                double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM],
                                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
    return -Zl*u*v;
}
