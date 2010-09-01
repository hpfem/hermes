.. _example-sing-pert:

Singular Perturbation
=====================

**Git reference:** Examples `singular perturbation <http://git.hpfem.org/hermes3d.git/tree/HEAD:/examples/singpert-aniso>`_.

We solve a singularly perturbed elliptic problem that exibits a thin anis-tropic boundary layer 
that is difficult to solve. This examples demonstrates how the anisotropic refinements can save 
a big amount of degrees of freedom.

.. index::
   single: mesh; dynamical
   single: problem; elliptic

The computational domain is the unit cube $(0, 1)^3$, and the equation solved has the form :

.. math::
   :nowrap:
   :label: sing-perturb

   \begin{eqnarray*}
   - \Delta u + K^2 u &= K^2 &\hbox{ in }\Omega \\
                    u &= 0 &\hbox{ on }\partial\Omega
   \end{eqnarray*}

The boundary conditions are homegeneous Dirichlet. The right-hand side is chosen to keep the 
solution $u(x,y,z) \approx 1$ inside the domain. Here, we choose $K^2 = 10^4$ but everything 
works for larger values of $K$ as well. 

It is quite important to perform the initial refinements towards the boundary, thus providing 
a better initial mesh for adaptivity and  making convergence faster. 

Convergence graphs:

.. image:: singpert-aniso-conv.png

.. image:: singpert-aniso-conv-time.png


Solution and hp-mesh:

.. image:: singpert-aniso-sln.png

.. image:: singpert-aniso-order.png

.. seealso::
  
   :ref:`example-elasto-statics`
