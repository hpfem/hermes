#include "definitions.h"

template <typename Scalar>
volume_matrix_acoustic_transient_planar_linear_form_1_1<Scalar>::volume_matrix_acoustic_transient_planar_linear_form_1_1(unsigned int i, unsigned int j, double ac_rho, double ac_vel)
    : MatrixFormVol<Scalar>(i, j), ac_rho(ac_rho), ac_vel(ac_vel)
{       
}

template <typename Scalar>
Scalar volume_matrix_acoustic_transient_planar_linear_form_1_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
                                          Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * (u->dx[i]*v->dx[i]+u->dy[i]*v->dy[i]);
    }
    return result / this->ac_rho;
}

template <typename Scalar>
Hermes::Ord volume_matrix_acoustic_transient_planar_linear_form_1_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
                                             Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
       result += wt[i] * (u->dx[i]*v->dx[i]+u->dy[i]*v->dy[i]);
    }	
    return result;
}

template <typename Scalar>
MatrixFormVol<Scalar>* volume_matrix_acoustic_transient_planar_linear_form_1_1<Scalar>::clone() const
{
    return new volume_matrix_acoustic_transient_planar_linear_form_1_1(*this);
}

template <typename Scalar>
volume_matrix_acoustic_transient_planar_linear_form_1_2<Scalar>::volume_matrix_acoustic_transient_planar_linear_form_1_2(unsigned int i, unsigned int j, double ac_rho, double ac_vel)
    : MatrixFormVol<Scalar>(i, j), ac_rho(ac_rho), ac_vel(ac_vel)
{       
}


template <typename Scalar>
Scalar volume_matrix_acoustic_transient_planar_linear_form_1_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
                                          Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += wt[i] *u->val[i] * v->val[i];
    }
    return result / (ac_rho * ac_vel *ac_vel * this->wf->get_current_time_step());
}

template <typename Scalar>
Hermes::Ord volume_matrix_acoustic_transient_planar_linear_form_1_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
                                             Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * u->val[i] * v->val[i];
    }	
    return result;
}

template <typename Scalar>
MatrixFormVol<Scalar>* volume_matrix_acoustic_transient_planar_linear_form_1_2<Scalar>::clone() const
{
    return new volume_matrix_acoustic_transient_planar_linear_form_1_2(*this);
}

template <typename Scalar>
volume_matrix_acoustic_transient_planar_linear_form_2_2<Scalar>::volume_matrix_acoustic_transient_planar_linear_form_2_2(unsigned int i, unsigned int j, double ac_rho, double ac_vel)
    : MatrixFormVol<Scalar>(i, j), ac_rho(ac_rho), ac_vel(ac_vel)
{       
}


template <typename Scalar>
Scalar volume_matrix_acoustic_transient_planar_linear_form_2_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
                                          Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * u->val[i] * v->val[i];
    }
    return -result;
}

template <typename Scalar>
Hermes::Ord volume_matrix_acoustic_transient_planar_linear_form_2_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
                                             Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
       result += wt[i] * u->val[i] * v->val[i];
    }	
    return -result;
}

template <typename Scalar>
MatrixFormVol<Scalar>* volume_matrix_acoustic_transient_planar_linear_form_2_2<Scalar>::clone() const
{
    return new volume_matrix_acoustic_transient_planar_linear_form_2_2(*this);
}

template <typename Scalar>
volume_matrix_acoustic_transient_planar_linear_form_2_1<Scalar>::volume_matrix_acoustic_transient_planar_linear_form_2_1(unsigned int i, unsigned int j, double ac_rho, double ac_vel)
    : MatrixFormVol<Scalar>(i, j), ac_rho(ac_rho), ac_vel(ac_vel)
{       
}

template <typename Scalar>
Scalar volume_matrix_acoustic_transient_planar_linear_form_2_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u,
                                          Hermes::Hermes2D::Func<double> *v, Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * u->val[i] * v->val[i];
    }
    return result / this->wf->get_current_time_step();
}

template <typename Scalar>
Hermes::Ord volume_matrix_acoustic_transient_planar_linear_form_2_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u,
                                             Hermes::Hermes2D::Func<Hermes::Ord> *v, Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
       result += wt[i] * u->val[i] * v->val[i];
    }	
    return result;
}

template <typename Scalar>
MatrixFormVol<Scalar>* volume_matrix_acoustic_transient_planar_linear_form_2_1<Scalar>::clone() const
{
    return new volume_matrix_acoustic_transient_planar_linear_form_2_1(*this);
}

template <typename Scalar>
volume_vector_acoustic_transient_planar_linear_form_1_2<Scalar>::volume_vector_acoustic_transient_planar_linear_form_1_2(unsigned int i, double ac_rho, double ac_vel)
    : VectorFormVol<Scalar>(i), ac_rho(ac_rho), ac_vel(ac_vel)
{
}

template <typename Scalar>
Scalar volume_vector_acoustic_transient_planar_linear_form_1_2<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
                                          Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;

    for (int i = 0; i < n; i++)
    {
        result += wt[i] * v->val[i] * ext[1]->val[i];
    }
    return result / (ac_rho * ac_vel * ac_vel * this->wf->get_current_time_step());
}

template <typename Scalar>
Hermes::Ord volume_vector_acoustic_transient_planar_linear_form_1_2<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                             Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * v->val[i] * ext[1]->val[i];
    }	
    return result;
}

template <typename Scalar>
VectorFormVol<Scalar>* volume_vector_acoustic_transient_planar_linear_form_1_2<Scalar>::clone() const
{
  return new volume_vector_acoustic_transient_planar_linear_form_1_2(*this);
}

template <typename Scalar>
volume_vector_acoustic_transient_planar_linear_form_2_1<Scalar>::volume_vector_acoustic_transient_planar_linear_form_2_1(unsigned int i, double ac_rho, double ac_vel)
    : VectorFormVol<Scalar>(i), ac_rho(ac_rho), ac_vel(ac_vel)
{
}

template <typename Scalar>
Scalar volume_vector_acoustic_transient_planar_linear_form_2_1<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *v,
                                          Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;

    for (int i = 0; i < n; i++)
    {
        result += wt[i] * v->val[i] * ext[0]->val[i];
    }
    return result / this->wf->get_current_time_step();
}

template <typename Scalar>
Hermes::Ord volume_vector_acoustic_transient_planar_linear_form_2_1<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                             Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
       result += wt[i] * v->val[i] * ext[0]->val[i];
    }	
    return result;
}

template <typename Scalar>
VectorFormVol<Scalar>* volume_vector_acoustic_transient_planar_linear_form_2_1<Scalar>::clone() const
{
  return new volume_vector_acoustic_transient_planar_linear_form_2_1(*this);
}

template <typename Scalar>
surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<Scalar>::surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance(unsigned int i, unsigned int j)
    : MatrixFormSurf<Scalar>(i, j)
{
}

template <typename Scalar>
Scalar surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<Scalar>::value(int n, double *wt, Hermes::Hermes2D::Func<Scalar> *u_ext[], Hermes::Hermes2D::Func<double> *u, Hermes::Hermes2D::Func<double> *v,
                                           Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::Func<Scalar> **ext) const
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * u->val[i] * v->val[i];
    }
    return result / ac_Z0;
}

template <typename Scalar>
Hermes::Ord surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<Scalar>::ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *u, Hermes::Hermes2D::Func<Hermes::Ord> *v,
                                              Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::Func<Hermes::Ord> **ext) const
{
    Hermes::Ord result(0);    
    for (int i = 0; i < n; i++)
    {
        result += wt[i] * u->val[i] * v->val[i];
    }	
    return result;

}

template <typename Scalar>
MatrixFormSurf<Scalar>* surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<Scalar>::clone() const
{
    return new surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance(*this);
}

template class volume_matrix_acoustic_transient_planar_linear_form_1_1<double>;
template class volume_matrix_acoustic_transient_planar_linear_form_1_2<double>;
template class volume_matrix_acoustic_transient_planar_linear_form_2_2<double>;
template class volume_matrix_acoustic_transient_planar_linear_form_2_1<double>;
template class volume_vector_acoustic_transient_planar_linear_form_1_2<double>;
template class volume_vector_acoustic_transient_planar_linear_form_2_1<double>;
template class surface_matrix_acoustic_transient_planar_linear_form_1_2_acoustic_impedance<double>;