NIST-12 (Multiple Difficulties)
-------------------------------

**Git reference:** Benchmark `nist-12 <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks/nist-12>`_.


This problem combines four difficulties of different strengths into the same problem by combining 
some of the features of the other test problems (reentrant corner, sharp peak, etc).


Model problem
~~~~~~~~~~~~~

Equation solved: Poisson equation

.. math::
    :label: NIST 12

       -\Delta u = f.

Domain of interest: L-shaped domain $(-1,1) \times (-1,1)$ \\ $(0,1) \times (-1,0)$.

Boundary conditions: Dirichlet, given by exact solution.

Right-hand side
~~~~~~~~~~~~~~~

Quite complicated, see below::

    {
    public:
      CustomRightHandSide(double alpha_p, double x_p, double y_p, double alpha_w, double x_w,
                          double y_w, double omega_c, double r_0, double epsilon)
        : DefaultNonConstRightHandSide(),
          alpha_p(alpha_p), x_p(x_p), y_p(y_p), alpha_w(alpha_w), x_w(x_w), y_w(y_w),
          omega_c(omega_c), r_0(r_0), epsilon(epsilon) { }

      double value(double x, double y) const {
        //For more elegant form please execute file "generate_rhs.py" 

        double a_P = (-alpha_p * pow((x - x_p), 2) - alpha_p * pow((y - y_p), 2));

        double a_W = pow(x - x_w, 2);
        double b_W = pow(y - y_w, 2);
        double c_W = sqrt(a_W + b_W);
        double d_W = ((alpha_w * x - (alpha_w * x_w)) * (2 * x - (2 * x_w)));
        double e_W = ((alpha_w * y - (alpha_w * y_w)) * (2 * y - (2 * y_w)));
        double f_W = (pow(alpha_w * c_W - (alpha_w * r_0), 2) + 1.0);
        double g_W = (alpha_w * c_W - (alpha_w * r_0));

        return -(4 * exp(a_P) * alpha_p * (alpha_p * (x - x_p) * (x - x_p) + alpha_p * (y - y_p) * (y - y_p) - 1)
               + ((alpha_w/(c_W * f_W)) - (d_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * d_W * g_W)/((a_W + b_W) * pow(f_W, 2)))
               + (alpha_w/(c_W * f_W)) - (e_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * e_W * g_W)/((a_W + b_W) * pow(f_W, 2))))
               + (1.0 / epsilon) * (1.0 / epsilon) * exp(-(1 + y) / epsilon));
      }


Exact solution
~~~~~~~~~~~~~~

.. math::

    u(x,y) =  r^{\alpha_{C} }\sin(\alpha_{C} \theta)
              + e^{-\alpha_{P} ((x - x_{P})^{2} + (y - y_{P})^{2})}
              + tan^{-1}(\alpha_{W} (r_{W} - r_{0}))  
              + e^{-(1 - y) / \epsilon}.

where $\alpha_C = \pi / \omega_C$, $r = \sqrt{x^2+y^2}$ and $\theta = tan^{-1}(y/x)$.  Here $\omega_C$ determines
the angle of the re-entrant corner. \

$(x_{P}, y_{P})$ is the location of the peak, $\alpha$ determines the strength of the peak, \

and $r_{W} = \sqrt{(x - x_{W})^{2} + (y - y_{W})^{2}}$. Here $(x_{W}, y_{W})$ is the center of the circular wave front.
$r_{0}$ is the distance from the wave front to the center of the circle, and $\alpha_W$ gives the steepness of the wave front. \

Last but not least, $\epsilon$ determines the strength of the boundary layer; the boundary layer was placed at $y = -1$.

Sample solution
~~~~~~~~~~~~~~~

Solution for $\omega_C = 3 \pi /2$,  $(x_{W}, y_{W}) = (0, -3/4)$,  $r_{0} = 3/4$, 
$\alpha_{W} = 200$,  $(x_{P}, y_{P}) = (\sqrt{5} / 4, -1/4)$,  $\epsilon = 1/100$:

.. image:: nist-12/solution.png
   :align: center
   :width: 600
   :alt: Solution.

Comparison of h-FEM (p=1), h-FEM (p=2) and hp-FEM with anisotropic refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (h-FEM, p=1, anisotropic refinements):

.. image:: nist-12/mesh_h1_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (h-FEM, p=2, anisotropic refinements):

.. image:: nist-12/mesh_h2_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-12/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-12/conv_dof_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-12/conv_cpu_aniso.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

hp-FEM with h-aniso and hp-aniso refinements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Final mesh (hp-FEM, h-anisotropic refinements):

.. image:: nist-12/mesh_hp_anisoh.png
   :align: center
   :width: 450
   :alt: Final mesh.

Final mesh (hp-FEM, hp-anisotropic refinements):

.. image:: nist-12/mesh_hp_aniso.png
   :align: center
   :width: 450
   :alt: Final mesh.

DOF convergence graphs:

.. image:: nist-12/conv_dof_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph.

CPU convergence graphs:

.. image:: nist-12/conv_cpu_hp.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

