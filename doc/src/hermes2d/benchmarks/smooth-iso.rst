Smooth-iso (Elliptic)
---------------------

**Git reference:** Benchmark `smooth-iso <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks/smooth-iso>`_.

We show that it is a very bad idea to approximate smooth solutions using low-order 
elements.

Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation 

.. math::
    :label: smooth-iso

       -\Delta u = f.

Domain of interest: Square $(0, \pi)^2$.

Right-hand side:

.. math::
    :label: smooth-iso-rhs
 
    f(x, y) = 2\sin(x)\sin(y).

Boundary conditions: Zero Dirichlet. 

In the code::

    // Enter boundary markers.
    BCTypes bc_types;
    bc_types.add_bc_dirichlet(BDY_DIRICHLET);

    // Enter Dirichlet boudnary values.
    BCValues bc_values;

Exact solution
~~~~~~~~~~~~~~

.. math::
    :label: smooth-iso-exact

    u(x, y) = \sin(x)\sin(y).

In the code::

    // Exact solution.
    static double fn(double x, double y)
    {
      return sin(x)*sin(y);
    }

    static double fndd(double x, double y, double& dx, double& dy)
    {
      dx = cos(x)*sin(y);
      dy = sin(x)*cos(y);
      return fn(x, y);
    }

Weak forms
~~~~~~~~~~

::

    // Weak forms.
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real>
    Real rhs(Real x, Real y)
    {
      return 2*sin(x)*sin(y);
    }

    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }

Sample solution
~~~~~~~~~~~~~~~

.. image:: benchmark-smooth-iso/sol_3d_view.png
   :align: center
   :width: 500
   :height: 300
   :alt: Solution.

Convergence comparisons
~~~~~~~~~~~~~~~~~~~~~~~

Below we show meshes obtained using various types of adaptivity. 
Note the tremendous differences in their performance. The meshes do not correspond to 
the same level of accuracy since the low-order methods could not achieve the same error 
as hp-FEM. Therefore, compare not only the number of DOF but also the error level. 
Convergence graphs for all cases are shown at the end of this section.

Final mesh (h-FEM, p=1): 27469 DOF, error 0.39173795799476 %

.. image:: benchmark-smooth-iso/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (h-FEM, p=2): 39185 DOF, error 0.0022127484879974 %

.. image:: benchmark-smooth-iso/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

Final mesh (hp-FEM): 49 DOF, error 4.2775412425017e-05 %

.. image:: benchmark-smooth-iso/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh

DOF convergence graphs:

.. image:: benchmark-smooth-iso/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: benchmark-smooth-iso/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.
