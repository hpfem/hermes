Heat Conduction
===============

**Git reference:** Examples `heat conduction <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes3d/examples/heat-conduction>`_.

This example describes the implementation of a simple time-dependent heat conduction problem inside the domain. 

.. index::
   single: mesh; fixed
   single: problem; time-depenedent
   single: problem; elliptic

The standard heat transfer PDE solved:

.. math::
   :nowrap:
   :label: heat-conduction

   \begin{eqnarray*}
   \frac{\partial T}{\partial t} - \Delta T &= & f \hbox{ in }\Omega \\ 
                                          T &= & 0 \hbox{ on }\partial\Omega
   \end{eqnarray*}

Domain of interest: Unit cube $(-1, 1)^3$:

.. image:: heat-conduction/heat-cond-domain.png

Right-hand side (load function):

.. math::
   :nowrap:
   :label: heat-conduction-rhs

   \begin{eqnarray*}
   f(x, y, z, t) & = & \mbox{cos}(t)(1-x^2)(1-y^2)(1-z^2)\ +\ \left(
                   2\mbox{sin}(t)(1-y^2)(1-z^2)
                   +2\mbox{sin}(t)(1-x^2)(1-z^2)
                   +2\mbox{sin}(t)(1-x^2)(1-y^2)
                   \right)
   \end{eqnarray*}

Equation :eq:`heat-conduction` is also equipped with an initial condition of the form: 

.. math::
   :label: heat-conduction-IC

   T(x, y, z, 0)  =  T_{init}(x, y, z) = 0  \ \ \ \mbox{in} \ \Omega. 

Exact solution is:

.. math:: 
   :nowrap:
   :label: heat-conduction-exact

   \begin{eqnarray*}
   T(x, y, z)  = \sin(t) (1 - x^2) (1 - y^2) (1 - z^2)
   \end{eqnarray*}

For simplicity we will use the implicit Euler method with a constant time step $\tau$, 
which transforms equation :eq:`heat-conduction` into: 

.. math::
   :label: heat-conduction-implicit

    \frac{T^{n+1} - T^n}{\tau} - \Delta T^{n+1} = 0.

The corresponding weak formulation is: 

.. math::
   :label: heat-conduction-form

    \int_{\Omega} \nabla T^{n+1}\cdot \nabla v + \int_{\Omega} \frac{T^{n+1}}{\tau} = 
    \int_{\Omega} f(t^{n+1}) v + \int_{\Omega} \frac{T^{n}}{\tau}.  

Code for the exact solution and the weak forms::

    double fn(double x, double y, double z)
    {
      return sin(TIME) * (1 - x*x) * (1 - y*y) * (1 - z*z);
    }

    double fndd(double x, double y, double z, double &dx, double &dy, double &dz)
    {
      dx = -2 * sin(TIME) * x * (1 - y*y) * (1 - z*z);
      dy = -2 * sin(TIME) * (1 - x*x) * y * (1 - z*z);
      dz = -2 * sin(TIME) * (1 - x*x) * (1 - y*y) * z;

      return fn(x, y, z);
     }

     // Boundary condition types.
     BCType bc_types(int marker) {
       return BC_ESSENTIAL;
     }

     template<typename real, typename scalar>
     scalar bilinear_form(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
     {
       return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e) + int_u_v<real, scalar>(n, wt, u, v, e) / TAU;
     }

     template<typename real> real rhs(real x, real y, real z)
     {
       real ddxx = -2 * sin(TIME) * (1 - y*y) * (1 - z*z);
       real ddyy = -2 * sin(TIME) * (1 - x*x) * (1 - z*z);
       real ddzz = -2 * sin(TIME) * (1 - x*x) * (1 - y*y);
       real dt = cos(TIME) * (1 - x*x) * (1 - y*y) * (1 - z*z);

       return dt - (ddxx + ddyy + ddzz);
     }

     template<typename real, typename scalar>
     scalar linear_form(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
     {
       return int_F_v<real, scalar>(n, wt, rhs, v, e) + int_u_v<real, scalar>(n, wt, data->ext + 0, v, e) / TAU;
     }

Before entering the main iteration loop, we need to initialize the previous solution sln_prev with the 
initial condition $T_{init}$ The solution class can be forced to return zero, to return a constant, 
or to return an arbitrary function using the methods set_zero(), set_const() and 
set_exact(), repectively. In this example, we initilize the temperature as all zero::

   // Construct initial solution and set zero.
   Solution sln_prev(&mesh);
   sln_prev.set_zero();

Next, the weak forms above are registered as following::

   // Initialize the weak formulation.

   WeakForm wf;
   wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
   wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY, &sln_prev);

Since the stiffness matrix does not depend on the solution, assembling only needs to be done once 
in the first time step. For all remaining time steps it will be the same, and we just need to 
re-construct the load vector. The code needs to be implemented. 

Solution graph:

.. image:: heat-conduction/heat-cond-sln.png
