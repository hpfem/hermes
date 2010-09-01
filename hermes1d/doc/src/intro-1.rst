============
Introduction
============

About Hermes1D
--------------

Hermes1D is an experimental C++ library for the solution of ordinary differential equations 
(ODE) and one-dimensional partial differential equations (PDE) with higher-order finite 
element methods (hp-FEM). In contrast to traditional time-stepping ODE solvers, Hermes1D 
constructs the solution using a variational principle. It starts from a weak formulation of 
the ODE/PDE problem and allows the equations to be defined in a very general implicit 
(vector-valued) form F(y, y', t) = 0. The approximation is a continuous, piecewise-polynomial 
function defined in the entire interval (0, T). In contrast to time-stepping schemes, the 
finite element approach makes it possible to prescribe boundary conditions either at the 
beginning or at the end of the time interval (combinations are possible for systems). The 
hp-FEM discretization leads to a system of nonlinear algebraic equations that is solved 
via the Newton's method or JFNK. Hermes1D comes 
with a free interactive online lab powered by UNR HPC cluster. The library is distributed 
under the BSD license. 

About this Document
-------------------

Prior to reading this document, we recommend that you install Hermes using instructions on 
its `home page <http://hpfem.org/hermes1d/>`_, and subscribe to the `mailing list 
<http://groups.google.com/group/hermes1d/>`_. Our mailing list is a very active place where 
you should get all answers quickly. 

The best way of reading this tutorial is to run the code at the same time. 
After making your way through the tutorial, you may want to browse the directory
with `examples <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/examples>`_ 
that contain a variety of different ODE and one-dimensional PDE  models. If you 
create an interesting model using Hermes, let us know and we will add it to the 
repository. 

The source code can be 
viewed in the `git repository <http://hpfem.org/git/gitweb.cgi/hermes1d.git/tree>`_.
For the 2D and 3D codes, see the `Hermes2D <http://hpfem.org/hermes2d/>`_ and 
`Hermes3D <http://hpfem.org/hermes3d/>`_ home pages, respectively.

User and Developer Documentation
--------------------------------

User documentation can be found in
the directory 'doc/'. Type 'make html' there to build it. The documentation is
available online at http://hpfem.org/hermes1d/doc/index.html.

To compile the C++ reference manual, go to 'hermes1d/doc.cpp/'. There
type 'doxygen hermes1d.lib-real.doxyfile'. The html files are in 
'h1d-real/html/index.html'. This documentation is also 
available online at http://hpfem.org/hermes1d/doc.cpp/h1d-real/html/index.html.

Mathematical Background
----------------------

When one speaks about the numerical solution of ODEs, one usually has in mind
initial value problems for equations of the form


.. math::

     {\d u_1\over\d x}=g_1(u_1, u_2, \dots, u_m, x),


.. math::
    :label: one

      \vdots


.. math::

     {\d u_m\over\d x}=g_m(u_1, u_2, \dots, u_m, x).

These are solved in a finite time interval $(0,T)$ using various time-stepping
methods. There are tons of those and some are quite sophisticated (meaning
multistep, higher-order, adaptive, etc.). But all of them have the following
common shortcomings:

* We would like to prescribe the initial value at $t = 0$ for some solution components and the end-time values at $t = T$ for others. Standard time stepping methods do not allow this.
* Global error control is problematic. One only can regulate the time step size locally -- this is something like "forward mesh refinement''. But one cannot do "backward mesh refinement'' or coarsening easily.
* We would like to prescribe a tolerance for the global error and then have the problem solved adaptively until this error tolerance is reached, without underresolving or overresolving too much. This is virtually impossible with adaptive time stepping methods.
* Standard time integration methods cannot change their order during the computation. For example, an adaptive RK4 method remains 4-order all the time. This is an analogy for $h$-refinement in FEM, and obviously it is highly inefficient. Correctly, the method should either do small low-order steps or large high-order steps to be efficient. We would like to see such an analogy of $hp$-refinement in ODE methods.
* We would like to solve more general ODEs than :eq:`one`.

This is why we decided to apply the $hp$-FEM methodology to ODEs and see what happens.

Equations
~~~~~~~~~

We implemented the first version of Hermes1D during one day while returning
from the 2009 SIAM CSE conference. First we considered the form :eq:`one` but
then we realized that with no extra work we can actually assume a much more
general implicit form


.. math::

     f_1(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x) = 0,


.. math::
    :label: two

      \vdots


.. math::

     f_m(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x) = 0.

Note that :eq:`two` contains :eq:`one` as a special case.
In fact, :eq:`two` can be written shortly as

.. math::
    :label: qqq

      \bfF(\bfU, \bfU', x) = 0

where ${\bfU} = (u_1,\dots,u_m)$ and ${\bfF} = (f_1,\dots,f_m)$.

Boundary conditions
~~~~~~~~~~~~~~~~~~~


So far, we have considered Dirichlet boundary conditions only, which can be
imposed either at the initial time $t = 0$ or the end-time $t = T$. Exactly one
condition per solution component has to be defined.


hp-FEM discretization
~~~~~~~~~~~~~~~~~~~~~


As always, the finite element discretization starts from a weak formulation.
With :eq:`two`, the situation is easy and we have


.. math::

     R_1(\bfY) = \int_0^T f_1(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_1 \, \d t = 0,


.. math::
    :label: three

      \vdots


.. math::

     R_N(\bfY) = \int_0^T f_m(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_N \, \d t = 0.

Here $v_1, v_2, \ldots, v_N$ are all basis functions for all solution
components (we can describe this more accurately if needed).  In the standard
sense, all basis functions corresponding to the solution component $u_i$ are
zero where $u_i$ has a Dirichlet boundary condition.  The vector $\bfY = (y_1,
y_2, \ldots, y_N)$ comprises all unknown coefficients of the finite element
basis functions for all solution components. The meshes for the solution
components $u_1, u_2, \ldots, u_m$ could (more precisely: *should*) be
different but for now we assume that they are the same.

Newton's method
~~~~~~~~~~~~~~~


We will drive the residual vector $\bfR = (R_1, R_2, \ldots, R_N)$ to zero
using the Newton's method. For that, we need the Jacobi matrix
$D\bfR/D\bfY$.

Let $1 \le i, j \le N$.
It is easy to calculate that

.. math::

     \frac{\partial R_i}{\partial y_j} = \int_0^T \frac{\partial f_{m(i)}}{\partial u_{n(j)}}(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_jv_i


.. math::
    :label: newt1

      + \frac{\partial f_{m(i)}}{\partial u'_{n(j)}}(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v'_jv_i \, \d t = 0.

Here, the function $m(i)$ takes a global index $1 \le i \le N$ and returns the
index of the function $f_{m(i)}$ which is associated with $R_i$. Analogously,
$n(j)$ takes a global index $1 \le j \le N$ and returns the index of the
solution component $u_{n(i)}$ where the basis function $v_j$ belongs to.

The integral in :eq:`newt1` has two parts because the functions $u_s$ and
$u'_s$ depend on the same solution coefficients.  Do not be confused by the
derivatives with respect to $u'_{n(j)}$ in :eq:`newt1`.  The functions $u_s$
and $u'_s$ are used as independent variables for the differentiation.


Interactive Web Accessibility
-----------------------------

* **Interactive web usage**. You can use Hermes (and other major open source FEM codes) remotely via any web browser, using the `FEMhub Online Numerical Methods Laboratory <http://lab.femhub.org/>`_. Your hardware will not be used as the online lab is powered by the University of Nevada, Reno (UNR) high-performance computing facility (`Research Grid <http://hpc.unr.edu/wiki/index.php/Main_Page>`_). You can compute with Hermes using an iPhone if you like.

.. image:: img/intro/iphone_large.png
   :align: center
   :width: 250
   :height: 450
   :alt: Hermes in iPhone.

See the `Hermes home page <http://hpfem.org/hermes1d/>`_ for more information. An overview of books, journal articles, conference proceedings papers and talks about Hermes and adaptive *hp*-FEM can be found in its `publications section <http://hpfem.org/publications/>`_.
