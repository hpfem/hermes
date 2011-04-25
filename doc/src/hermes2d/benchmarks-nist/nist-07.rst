NIST-07 (Boundary Line Singularity)
-----------------------------------

**Git reference:** Benchmark `nist-07 <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks/nist-07>`_.

Many papers on testing adaptive algorithms use a 1D example with a singularity of the form $x^{\alpha}$
at the left endpoint of the domain. This can be extended to 2D by simply making the solution be
constant in $y$.


Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation 

.. math::
    :label: line-sing

       -\Delta u = f.

Domain of interest: Unit Square $(0, 1)^2$

Boundary conditions: Dirichlet, given by exact solution.

Exact solution
~~~~~~~~~~~~~~

.. math::

    u(x,y) = x^{\alpha} 

where $\alpha \geq 0.5$ determines the strength of the singularity.

Right-hand side: Obtained by inserting the exact solution into the equation.
The corresponding code snippet is shown below::

    {
    public:
      CustomRightHandSide(double alpha) : DefaultNonConstRightHandSide(), alpha(alpha) {};

      virtual double value(double x, double y) const {
        return - alpha * (alpha - 1.) * pow(x, alpha - 2.);
    } 


Sample solution
~~~~~~~~~~~~~~~

Solution for $\alpha = 0.6$:

.. image:: nist-07/solution.png
   :align: center
   :width: 600
   :height: 400
   :alt: Solution.

Comparison of h-FEM (p=1), h-FEM (p=2) and hp-FEM with anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM, p=1, anisotropic refinements):

.. image:: nist-07/mesh_h1_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (h-FEM, p=2, anisotropic refinements):

.. image:: nist-07/mesh_h2_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-07/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-07/conv_dof_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-07/conv_cpu_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

hp-FEM with h-aniso and hp-aniso refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-07/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, hp-anisotropic refinements):

.. image:: nist-07/mesh_hp_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-07/conv_dof_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-07/conv_cpu_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

