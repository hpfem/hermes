Interior Layer (Elliptic)
-------------------------

**Git reference:** Benchmark `layer-interior <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks/layer-interior>`_.

This example has a smooth solution that exhibits a steep interior layer.

Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation 

.. math::
    :label: layer-interior

       -\Delta u = f.

Domain of interest: Unit square $(0, 1)^2$.

Right-hand side:

.. math::
    :label: layer-interior-rhs
 
    f(x, y) = \frac{27}{2} (2y + 0.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2} +
    \frac{27}{2} (2x - 2.5)^2 (\pi - 3t) \frac{S^3}{u^2 t_2}
    - \frac{9}{4} (2y + 0.5)^2 \frac{S}{u t^3} -
    \frac{9}{4} (2x - 2.5)^2 \frac{S}{u t^3} +
    18 \frac{S}{ut}.

Exact solution
~~~~~~~~~~~~~~

.. math::
    :label: layer-interior-exact

    u(x, y) = \mbox{atan}\left(S \sqrt{(x-1.25)^2 + (y+0.25)^2} - \pi/3\right).

where $S$ is a parameter (slope of the layer). With larger $S$, this problem 
becomes difficult for adaptive algorithms, and at the same time the advantage of 
adaptive $hp$-FEM over adaptive low-order FEM becomes more significant. We will 
use $S = 60$ in the following.

In the code::

    // Exact solution.
    static double fn(double x, double y)
    {
      return atan(SLOPE * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
    }
    
    static double fndd(double x, double y, double& dx, double& dy)
    {
      double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
      double u = t * (sqr(SLOPE) * sqr(t - M_PI/3) + 1);
      dx = SLOPE * (x-1.25) / u;
      dy = SLOPE * (y+0.25) / u;
      return fn(x, y);
    }
    
Boundary conditions
~~~~~~~~~~~~~~~~~~~

Nonconstant Dirichlet, dictated by the choice of the exact solution::

    // Essential (Dirichlet) boundary condition values.
    scalar essential_bc_values(double x, double y)
    {
      return fn(x, y);
    }

The callback is registered as follows::

    // Enter boundary markers.
    BCTypes bc_types;
    bc_types.add_bc_dirichlet(BDY_DIRICHLET);

    // Enter Dirichlet boudnary values.
    BCValues bc_values;
    bc_values.add_function(BDY_DIRICHLET, essential_bc_values);

Weak forms
~~~~~~~~~~
    
::

    // Bilinear form for the Poisson equation.
    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }
    
    template<typename Real>
    Real rhs(Real x, Real y)
    {
      Real t2 = sqr(y + 0.25) + sqr(x - 1.25);
      Real t = sqrt(t2);
      Real u = (sqr(M_PI - 3.0*t)*sqr(SLOPE) + 9.0);
      return 27.0/2.0 * sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) +
             27.0/2.0 * sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * pow(SLOPE,3.0) / (sqr(u) * t2) -
             9.0/4.0 * sqr(2.0*y + 0.5) * SLOPE / (u * pow(t,3.0)) -
             9.0/4.0 * sqr(2.0*x - 2.5) * SLOPE / (u * pow(t,3.0)) +
             18.0 * SLOPE / (u * t);
    }
     
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
    }

Sample solution
~~~~~~~~~~~~~~~

.. image:: benchmark-layer-interior/sol_3d_view.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Convergence comparison
~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM with linear elements):

.. image:: benchmark-layer-interior/mesh-h1.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with linear elements).

Final mesh (h-FEM with quadratic elements):

.. image:: benchmark-layer-interior/mesh-h2.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (h-FEM with quadratic elements).

Final mesh (hp-FEM):

.. image:: benchmark-layer-interior/mesh-hp.png
   :align: center
   :width: 500
   :height: 400
   :alt: Final mesh (hp-FEM).

DOF convergence graphs:

.. image:: benchmark-layer-interior/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU time convergence graphs:

.. image:: benchmark-layer-interior/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.
