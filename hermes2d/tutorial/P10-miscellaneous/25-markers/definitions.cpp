class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double A_SE, double A_NE, double A_SW, double A_NW, double rhs) : WeakForm(1) {
    add_matrix_form(new MatrixFormVol(0, 0, A_SE, A_NE, A_SW, A_NW));
    add_vector_form(new VectorFormVol(0, rhs));
    add_vector_form_surf(new VectorFormSurf(0, A_SE, A_NE, A_SW, A_NW));
  };

private:
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j, double A_SE, double A_NE, double A_SW, double A_NW) 
                  : WeakForm::MatrixFormVol(i, j), A_SE(A_SE), A_NE(A_NE), A_SW(A_SW), A_NW(A_NW) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) == SOUTH_EAST)
        return this->A_SE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
      if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) == NORTH_EAST)
        return this->A_NE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
      if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) == SOUTH_WEST)
        return this->A_SW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
      if(wf->get_element_markers_conversion()->get_user_marker(e->elem_marker) == NORTH_WEST)
        return this->A_NW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }

    double A_SE, A_NE, A_SW, A_NW;
  };

  class VectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVol(int i, double rhs) : WeakForm::VectorFormVol(i), rhs(rhs) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      return rhs * int_v<Real>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    double rhs;
  };

  class VectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurf(int i, double A_SE, double A_NE, double A_SW, double A_NW) 
                   : WeakForm::VectorFormSurf(i), A_SE(A_SE), A_NE(A_NE), A_SW(A_SW), A_NW(A_NW) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_VERTICAL_SE)
        return - A_SE * int_v<Real>(n, wt, v);
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_VERTICAL_NE)
        return - A_NE * int_v<Real>(n, wt, v);
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_VERTICAL_NW)
        return - A_NW * int_v<Real>(n, wt, v);
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_VERTICAL_SW)
        return - A_SW * int_v<Real>(n, wt, v);
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_TOP_NE)
        return A_NE * int_v<Real>(n, wt, v);
      if(wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker) == BDY_TOP_NW)
        return A_NW * int_v<Real>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return int_v<Ord>(n, wt, v);
    }

    double A_SE, A_NE, A_SW, A_NW;
  };

};
