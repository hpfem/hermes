General 2nd-Order Linear Equation (07-general)
--------------------------------------

**Git reference:** Tutorial example `07-general <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P01-linear/07-general>`_. 

Model problem
~~~~~~~~~~~~~

This example deals with a linear second-order equation of the form 

.. math::

         -\frac{\partial}{\partial x}\left(a_{11}(x,y)\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial x}\left(a_{12}(x,y)\frac{\partial u}{\partial y}\right) - \frac{\partial}{\partial y}\left(a_{21}(x,y)\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial y}\left(a_{22}(x,y)\frac{\partial u}{\partial y}\right) + a_1(x,y)\frac{\partial u}{\partial x} + a_{21}(x,y)\frac{\partial u}{\partial y} + a_0(x,y)u = rhs(x,y),

equipped with Dirichlet and/or Neumann boundary conditions. Its goal is to show how to 
use space-dependent coefficients and how to define quadrature orders explicitly. 

Defining non-constant equation coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we define the (generally) non-constant equation coefficients:
::

    double a_11(double x, double y) {
      if (y > 0) return 1 + x*x + y*y;
      else return 1;
    }

and so on. Then we define boundary conditions as usual. 

Manual setting of integration orders in weak forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weak formulation contains both volumetric and surface integrals. 

The Ord class in Hermes (see the file `src/weakform/forms.h 
<http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/weakform/forms.h>`_) provides
an automatic parser of weak forms that is able to determine the integration orders for 
algebraic expressions. So, in order to define an integration order explicitly, one can 
provide on top the weak form another function that defines a simple algebraic expression 
that leads the parser to the desired polynomial degree. The values defined in this  
additional function are not used for computation. 

::

    // (Volumetric) bilinear form
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i=0; i < n; i++) {
        double x = e->x[i];
        double y = e->y[i];
        result += (a_11(x, y)*u->dx[i]*v->dx[i] + 
                   a_12(x, y)*u->dy[i]*v->dx[i] +
                   a_21(x, y)*u->dx[i]*v->dy[i] +
                   a_22(x, y)*u->dy[i]*v->dy[i] +
                   a_1(x, y)*u->dx[i]*v->val[i] +
                   a_2(x, y)*u->dy[i]*v->val[i] +
                   a_0(x, y)*u->val[i]*v->val[i]) * wt[i];
      }
      return result;
    }

    // Integration order for the bilinear form
    Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis 
                                                        // and test function plus two
    }

    // Surface linear form (natural boundary conditions)
    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, g_N, v, e);
    }
  
    // Integration order for surface linear form
    Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }
  
    // Volumetric linear form (right-hand side)
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }
  
    // Integration order for the volumetric linear form
    Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
    }

Note the sign of the surface linear form. When solving a linear problems, i.e., when initializing the 
DiscreteProblem class with the third parameter being is_linear = true, all linear forms have to be on 
the right-hand side and all bilinear forms on the left. 

Sample results
~~~~~~~~~~~~~~

The output of this example is shown below:

.. image:: 07/general.png
   :align: center
   :width: 500
   :height: 400
   :alt: Output of example 07-general.
