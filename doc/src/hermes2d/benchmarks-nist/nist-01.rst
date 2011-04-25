NIST-01 (Analytic Solution)
---------------------------

**Git reference:** Benchmark `nist-01 <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks/nist-01>`_.

The purpose of this benchmark is to observe how an adaptive algorithm behaves in a context where 
adaptivity isn’t really needed (smooth solution). 


Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation 

.. math::
    :label: NIST-1

       -\Delta u = f.

Domain of interest: Unit Square $(0, 1)^2$.

Boundary conditions: Dirichlet, given by exact solution.

Exact solution
~~~~~~~~~~~~~~

.. math::

    u(x,y) = 2^{4p}x^{p}(1-x)^{p}y^{p}(1-y)^p

where parameter $p$ determines the degree of the polynomial solution. 

Right-hand side: Obtained by inserting the exact solution into the latter equation.
The corresponding code snippet is shown below::

    {
    public:
      CustomRightHandSide(double pol_deg)
        : DefaultNonConstRightHandSide(), pol_deg(pol_deg) {};

      virtual double value(double x, double y) const {
        double a = pow(2.0, 4.0*pol_deg);
        double b = pow(x-1.0, 8.0);
        double c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
        double d = pow(y-1.0, pol_deg);
        double e = pow(y-1.0, 8.0);
        double f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
        double g = pow(x-1.0, pol_deg);

        return -(pol_deg*a*pow(x, 8.0)*b*c*pow(y, pol_deg)*d
             + pol_deg*a*pow(y, 8.0)*e*f*pow(x, pol_deg)*g);
    }


Sample solution
~~~~~~~~~~~~~~~

Solution for $p = 10$:

.. image:: nist-01/solution.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Comparison of h-FEM (p=1), h-FEM (p=2) and hp-FEM with anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM, p=1, anisotropic refinements):

.. image:: nist-01/mesh_h1_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (h-FEM, p=2, anisotropic refinements):

.. image:: nist-01/mesh_h2_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-01/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-01/conv_dof_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-01/conv_cpu_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

hp-FEM with iso, h-aniso and hp-aniso refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (hp-FEM, isotropic refinements):

.. image:: nist-01/mesh_hp_iso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-01/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, hp-anisotropic refinements):

.. image:: nist-01/mesh_hp_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-01/conv_dof_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-01/conv_cpu_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.


