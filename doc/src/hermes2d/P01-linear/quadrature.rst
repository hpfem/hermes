Determination of Quadrature Orders in Weak Forms
------------------------------------------------

You may wonder why templates are used in the definition of weak forms. As a matter of fact, 
they do not have to be, as we will see in a moment. However, if the weak form only contains 
algebraic operations (without if-then statements and such), templates help to determine
numerical integration orders automatically. In higher-order FEM, basis and test functions may 
have very different polynomial degrees, ranging from one and some maximum polynomial 
degree (currently 10 in Hermes). The basis and test functions can be combined inside the 
weak forms in many different ways. As a result, the minimum quadrature order which is needed 
to evaluate a weak form accurately may vary between zero (product of gradients of 
two linear functions) to infinity (whenever a nonpolynomial expression is present). 
Numerical quadrature is one of the trickiest issues in higher-order FEM.

Brute force does not work
~~~~~~~~~~~~~~~~~~~~~~~~~

A brute-force solution to this problem would be to integrate everything using 
a maximum order, but this would lead to tremendous computing times. Therefore Hermes offers 
two options: the polynomial degree of the integrated expressions can be detected 
automatically (via templates), or the user can define for each weak form the 
quadrature order manually. If the weak form only contains polynomial expressions, 
the former approach works very well. If the form is more complicated, it is recommended 
to handle the integration orders manually. 

Automatic determination of quadrature order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In example 03-poisson, the bilinear and linear forms were defined using templates,

::

    // return the value \int \nabla u . \nabla v dx
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    // return the value \int v dx
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return CONST_F * int_v<Real, Scalar>(n, wt, v);
    }

and registered using the callback() macro,

::

    // initialize the weak formulation
    WeakForm wf();
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_vector_form(callback(linear_form));
   
The callback() macro, defined in `src/weakform/forms.h 
<http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/weakform/forms.h>`_ by

::

    #define callback(a)     a<double, scalar>, a<Ord, Ord>

expands the above add_matrix_form() and add_vector_form() functions into

::

    // initialize the weak formulation
    WeakForm wf();
    wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>);
    wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);

For those who are not familiar with templates, they make it possible to 
call the same function with different parameter types. In particular, 
using bilinear_form<double, scalar> and bilinear_form<Ord, Ord> for
the bilinear form defined above gives 

::

    scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return int_grad_u_grad_v<double, scalar>(n, wt, u, v);
    }

    Ord bilinear_form(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }

The <double, scalar> copy is used to obtain the result of the numerical integration,
the <Ord, Ord> copy for automatic evaluation of the quadrature order. 
The parser (see `src/weakform/forms.h 
<http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/src/weakform/forms.h>`_) 
works well for algebraic expressions. If the weak form bilinear_form() is complicated, 
one can create and register a simpler weak form bilinear_form_order() for the parser,
that provides an arbitrary expression with the same polynomial degree as 
the integrand in bilinear_form(). Then the two functions would be registered as 

::

    wf.add_matrix_form(bilinear_form, bilinear_form_order);

Of course the same holds for linear forms.
If the bilinear form contains things like the if-then statement, it cannot 
be parsed. Whenever the weak form contains non-polynomial expressions or
is otherwise very complicated, it is recommended to handle the quadrature 
orders manually.

Manual determination of quadrature order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The polynomial degree of basis and test functions inside a bilinear or linear form 
can be handled manually as follows

::

    Ord bilinear_form_order(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      int uo = u->val[0].get_order();
      int vo = v->val[0].get_order();
      return Ord(uo + vo);            // this would correspond to integral of u times v
    }

It is also possible to return a constant order (for example 5) by using 

::

    Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return Ord(5);
    }

Currently, one cannot make the integration order dependent on spatial coordinates and such. However,
one can assign different weak forms to elements with different material markers. This is
described in examples `neutronics-iron-water <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/neutronics-iron-water.html>`_,
`neutronics-saphir <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/examples/neutronics-saphir.html>`_ and others.

The following example 07-general handles quadrature orders manually. 
