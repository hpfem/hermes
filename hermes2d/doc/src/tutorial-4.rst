=========================================================
Tutorial Part IV (Adaptivity for Time-Dependent Problems)
=========================================================

(Space-time) adaptive FEM and *hp*-FEM for time-dependent PDE and PDE systems is one of 
the most advanced techniques Hermes can do. Although we have published it 
in several `scientific articles 
<http://hpfem.math.unr.edu/people/pavel/public/papers.html>`_, 
there is still a lot of space for improvement. Let us know through the
`Hermes2D mailing list <http://groups.google.com/group/hermes2d/>`_ if 
you are interested in getting involved. In this part of the tutorial 
we explain the basic idea of the method and show several examples.

Adaptive hp-FEM with Dynamical Meshes
-------------------------------------

The adaptivity with dynamical meshes in Hermes is based on the combination 
of the multimesh *hp*-FEM with the classical Rothe's method. 

The Rothe's method is a natural counterpart of the widely used Method of Lines (MOL). 
Recall that the MOL performs discretization in space while 
keeping the time variable continuous, which leads to a system of ODEs in time. The Rothe's 
method, on the contrary, preserves the continuity of the spatial variable while discretizing time. 
In every time step, an evolutionary PDE is approximated by means of one or more time-independent ones. 
The number of the time-independent equations per time step is proportional to the order of accuracy of the 
time discretization method. For example, when employing the implicit Euler method, one 
has to solve one time-independent PDE per time step. The Rothe's method is fully equivalent to the 
MOL if no adaptivity in space or time takes place, but it provides a better setting 
for the application of spatially adaptive algorithms. The spatial discretization error
can be controlled by solving the time-independent equations adaptively, and the size of 
the time step can be adjusted using standard ODE techniques. 

For the sake of clarity, let us consider a simple linear parabolic problem 

.. math::

    \frac{\partial u}{\partial t} - \Delta u = f

and discretize the time variable using the implicit Euler method. We obtain 

.. math::

    \frac{u^{n+1} - u^n}{\tau} - \Delta u^{n+1} = f^{n+1},

where 

.. math::

    u^{n+1}(x,y) \approx u(x, y, t^{n+1})\ \mbox{and} \  f^{n+1}(x, y) \approx f(x, y, t^{n+1}).

The equation for $u^{n+1}$ does no longer depend on time and we can solve it adaptively 
as any other time-independent equation (or equation system). The only thing worth 
mentioning is that the previous time step approximation $u^n$ is now defined on 
a locally refined mesh that was obtained during the previous time step. This 
situation, however, can be handled routinely via the multimesh discretization 
method. In fact, the user does not have to worry about anything. The methodology is 
illustrated below.

Nonlinear Parabolic Problem (22)
--------------------------------

**Git reference:** Tutorial example `22-newton-timedep-heat-adapt 
<http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/22-newton-timedep-heat-adapt>`_.

Description coming soon.


Heat and Moisture Transfer in Concrete (23)
-------------------------------------------

**Git reference:** Tutorial example `23-heat-and-moisture-adapt 
<http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/23-heat-and-moisture-adapt>`_.

Description coming soon.

Gross-Pitaevski Equation (24)
-----------------------------

**Git reference:** Tutorial example `24-newton-timedep-gp-adapt 
<http://git.hpfem.org/hermes2d.git/tree/HEAD:/tutorial/24-newton-timedep-gp-adapt>`_.

Description coming soon.
