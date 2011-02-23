Adaptive low-order FEM and hp-FEM
---------------------------------

In traditional low-order FEM, refining an element is not algorithmically complicated,
and so the most difficult part is to find out what elements should be
refined. To do this, people employ various techniques ranging from rigorous
guaranteed a-posteriori error estimates to heuristic criteria such as residual
error indicators, error indicators based on steep gradients, etc. However, 
these approaches are not suitable for multiphysics coupled problems nor for 
higher-order finite element methods: Rigorous guaranteed error
estimates only exist for very simple problems (such as linear elliptic PDE),
and moreover only for low-order finite elements. 
Note that virtually no a-posteriori error estimates capable of 
guiding automatic hp-adaptivity are available even for simplest elliptic problems. 
The heuristic techniques listed above are not employed in Hermes since they may fail 
in non-standard situations, and because they lack a transparent relation to the 
true approximation error.

Why is low-order FEM so inefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adaptive low-order FEM is known to be notoriously inefficient, and practitioners
are rightfully skeptical of it. The reason is its extremely slow convergence 
that makes large computations freeze. 
This is illustrated in the following graphs that compare typical convergences 
of adaptive FEM with linear elements, adaptive FEM with quadratic elements, and 
adaptive hp-FEM:

.. image:: conv-intro/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: Typical convergence curves for adaptive linear FEM, quadratic FEM, and *hp*-FEM.

The reader can see that the 
linear FEM would need in the order of 1,000,000,000,000,000,000 degrees of freedom 
(DOF) to reach a level of accuracy where the hp-FEM is with less than 10,000 DOF. 
A similar effect can be observed in the CPU-time convergence graph:

.. image:: conv-intro/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: Typical convergence curves for adaptive linear FEM, quadratic FEM, and *hp*-FEM.

These convergence curves are typical representative examples, confirmed with
many numerical experiments of independent researchers, and supported with
theory. The low-order FEM is doomed by the underlying math -- its poor convergence cannot 
be fixed by designing smarter adaptivity algorithms.

Why is hp-FEM so efficient
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to obtain fast, usable adaptivity (the red curve), one
has to resort to adaptive *hp*-FEM. The *hp*-FEM takes advantage of 
the following facts:

* Large high-degree elements approximate smooth parts of a solution *much* 
  better than small linear ones. 
  We have a benchmark `smooth-iso <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks/smooth-iso.html>`_ 
  that illustrates this in an impressive way - please spend a few minutes to check it out, 
  the results will surprize you. In the 
  Hermes2D repository, it can be found in the directory 
  `benchmarks/ <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/benchmarks>`_.
* This holds the other way round for singularities,
  steep gradients, oscillations and other "bad" spots. These are 
  approximated best using small low-order elements.
* In order to capture efficiently anisotropic solution behavior, one needs adaptivity algorithms 
  that can refine meshes anisotropically both in $h$ and $p$. This is illustrated 
  in  benchmarks 
  `smooth-aniso-x <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks/smooth-aniso-x.html>`_  
  and `nist-07(line singularity) <http://hpfem.org/hermes/doc/src/hermes2d/nist/nist-07.html>`_.

What it takes to do adaptive hp-FEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatic adaptivity in the *hp*-FEM is substantially different from adaptivity
in low-order FEM, since every element can be refined in many different ways.
The following figure shows several illustrative refinement candidates for 
a fourth-order element.

.. image:: conv-intro/refinements.png
   :align: center
   :width: 650
   :height: 300
   :alt: Examples of *hp*-refinements.

Of course, the number of possible element refinements is implementation-dependent.
In general it is very low in $h$ or $p$ adaptivity, much higher in $hp$ adaptivity, 
and it rises even more when anisotropic refinements are enabled. This is why Hermes 
has eight different adaptivity options P_ISO, P_ANISO, H_ISO, H_ANISO,
HP_ISO, HP_ANISO_P, HP_ANISO_H, HP_ANISO. In this order, usually P_ISO yields the 
worst results and HP_ANISO the best. However, even P_ISO can be very efficient with 
a good starting mesh.  In the most general HP_ANISO 
option, around 100 refinement candidates for each element are considered. 
Naturally, the adaptivity algorithm takes progressively more time as more 
refinement candidates are probed. The difference between the HP_ANISO_H
option (next best to HP_ANISO) and HP_ANISO is quite significant. So, this is 
where the user can use his a-priori knowledge of the solution to make the 
computation faster. 

Why do we need more than standard error estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to the large number of refinement options in each element, classical error estimators that
provide just one number per element are not enough. To guide hp-adaptivity, one really needs 
to know the **shape** of the approximation error, not only its magnitude.

In analogy to the most successful adaptive ODE solvers,
Hermes uses a pair of approximations with different orders of accuracy to obtain
this information: *coarse mesh solution* and 
*fine mesh solution*. The initial coarse mesh is read from the mesh file,
and the initial fine mesh is created through its global refinement both in
$h$ and $p$.

The fine mesh solution is the approximation of interest both during the adaptive
process and at the end of computation. After years of experimentation, the coarse 
mesh solution was gradually replaced in Hermes with global orthogonal projection of the fine 
mesh solution on the coarse mesh. In most cases, this yields a better 
convergence behavior than using the coarse mesh solve, and the projection 
problem is always linear and well conditioned. 

Robustness of the reference solution approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the reference solution approach is PDE independent, which is truly great for multiphysics
coupled problems. Hermes does not use a single analytical error estimate 
or any other technique that would narrow down its applicability to selected 
equations or low-order FEM. 

Room for improvement
~~~~~~~~~~~~~~~~~~~~

An obvious disadvantage of the reference solution approach to automatic adaptivity is its higher 
computational cost, especially in 3D. We are aware of this fact and would not mind 
at all replacing the current paradigm  
with some cheaper technique -- as long as it is PDE-independent, 
works for elements of high orders, and handles anisotropy in both 'h' and 'p'.
Seemingly, however, no such alternatives exist. 
