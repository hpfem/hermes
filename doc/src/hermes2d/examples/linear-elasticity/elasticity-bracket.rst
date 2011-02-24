Bracket (Linear Elasticity)
---------------------------

**Git reference:** Example `elasticity-bracket <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/linear-elasticity/elasticity-bracket>`_.

This example employs the adaptive multimesh hp-FEM to solve equations of linear elasticity. The code is
basically the same as in `example crack <file:///home/pavel/repos/hermes/doc/_build/html/src/hermes2d/examples/linear-elasticity/elasticity-crack.html>`_
and thus we do not discuss it in great detail.

Model problem
~~~~~~~~~~~~~

Our domain is a bracket loaded on its top edge and fixed to the wall:

.. math::
    :nowrap:

    \begin{eqnarray*}   \bfu \!&=&\! 0 \ \ \ \ \ \rm{on}\ \Gamma_1  \\   \dd{u_2}{n} \!&=&\! f \ \ \ \ \ \rm{on}\ \Gamma_2 \\   \dd{u_1}{n} = \dd{u_2}{n} \!&=&\! 0 \ \ \ \ \ \rm{elsewhere.} \end{eqnarray*}

The dimensions are L = 0.7 m, T = 0.1 m and the force $f = 10^3$ N.

.. image:: example-bracket/bracket.png
   :align: center
   :width: 400
   :height: 400
   :alt: Computational domain for the elastic bracket problem.

Sample results
~~~~~~~~~~~~~~

The following figures show the two meshes and their polynomial
degrees after several adaptive steps: 

.. image:: example-bracket/sys-xorders.png
   :align: left
   :width: 300
   :height: 300
   :alt: $x$ displacement -- mesh and polynomial degrees.

.. image:: example-bracket/sys-yorders.png
   :align: right
   :width: 300
   :height: 300
   :alt: $y$ displacement -- mesh and polynomial degrees.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">


Note that the meshes are slightly different, not only in
polynomial degrees, but also in element refinements. This is 
possible in Hermes thanks to 
`adaptive multi-mesh hp-FEM <file:///home/pavel/repos/hermes/doc/_build/html/src/hermes2d/linear-adapt/multimesh.html>`_.

Convergence comparison
~~~~~~~~~~~~~~~~~~~~~~

Convergence graphs of adaptive h-FEM with linear elements, h-FEM with quadratic elements
and hp-FEM are shown below.

.. image:: example-bracket/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for tutorial example 11-adapt-system.

The following graph shows convergence in terms of CPU time. 

.. image:: example-bracket/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for example bracket

Comparison of the multimesh and single-mesh hp-FEM: 

.. image:: example-bracket/conv_compar_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: comparison of multimesh and single mesh hp-FEM

.. image:: example-bracket/conv_compar_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: comparison of multimesh and single mesh hp-FEM

In this example the difference between the multimesh *hp*-FEM and the single-mesh
version was not extremely large since the two elasticity equations are very 
strongly coupled and have singularities at the same points. 
To see more significant differences, look at the tutorial 
example `P04-linear-adapt/02-system-adapt <file:///home/pavel/repos/hermes/doc/_build/html/src/hermes2d/linear-adapt/multimesh-example.html>`_.
