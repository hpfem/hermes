====================
Hermes Documentation
====================

Thank you for your interest in Hermes!

Hermes is a C++ library for rapid development of adaptive *hp*-FEM and *hp*-DG solvers.

This document is organized as follows: 

* Section 1 provides general information about Hermes and the computational methods it uses,
  and gives references to underlying scientific articles.
* Section 2 describes how to install Hermes on various hardware platforms, and how to 
  install matrix solver packages and various optional packages. 
* Section 3 explains how to use Git and Github, and how you can contribute to the project if interested.
* Section 4 contains a tutorial to Hermes2D. Please read this tutorial first even if you are 
  interested in Hermes3D or Hermes1D, since their usage is virtually the same. The tutorial 
  will walk you in small steps through the solution
  of linear, nonlinear, and time-dependent problems from various engineering and scientific areas, 
  using higher-order elements and adaptivity algorithms, and solving multiphysics coupled problems. 
* Section 5 shows how Hermes performs on numerous benchmarks with known exact solutions. Many 
  of them come from the National Institute for Standards and Technology (NIST). Benchmarks 
  like this are great for assessing the performance of a finite element code.
* Section 6 presents examples from various application areas such as acoustics, fluid and solid
  mechanics, electromagnetics, neutronics, quantum chemistry, ground-water flow, and others. 
* Section 7 describes the usage of Hermes3D. Since it is very similar to Hermes2D, 
  we mostly focus on their differences.
* Section 8 shows how to use Hermes1D. Again, we mainly point out where it differs
  from Hermes2D.

This document is under continuous development. If you find bugs, typos, dead links 
and such, please report them to one of the mailing lists for 
`Hermes1D <http://groups.google.com/group/hermes1d/>`_,
`Hermes2D <http://groups.google.com/group/hermes2d/>`_, or
`Hermes3D <http://groups.google.com/group/hermes3d/>`_ -- thanks!

1. Introduction
---------------

.. toctree::
    :maxdepth: 1

    src/about-hermes
    src/math-background
    src/web-access
    src/citing-hermes

2. Installation
---------------

.. toctree::
    :maxdepth: 1

    src/installation/linux
    src/installation/mac
    src/installation/win-cygwin
    src/installation/win-msvc
    src/installation/matrix_solvers
    src/installation/cython_installation
    src/installation/exodusII_netcdf

3. Collaboration
----------------

.. toctree::
    :maxdepth: 1

    src/first_pull_request

4. Tutorial
-----------

.. toctree::
    :maxdepth: 1

    src/hermes2d/P01-linear
    src/hermes2d/P02-nonlinear
    src/hermes2d/P03-timedep
    src/hermes2d/P04-linear-adapt
    src/hermes2d/P05-nonlinear-adapt
    src/hermes2d/P06-timedep-adapt
    src/hermes2d/P07-eigen
    src/hermes2d/P08-fvm-and-dg
    src/hermes2d/P09-trilinos
    src/hermes2d/P10-miscellaneous

5. Benchmarks
-------------

.. toctree::
    :maxdepth: 1

    src/hermes2d/benchmarks-nist
    src/hermes2d/benchmarks-general

6. Examples
-----------

.. toctree::
    :maxdepth: 1

    src/hermes2d/examples/acoustics.rst
    src/hermes2d/examples/advection-diffusion-reaction.rst
    src/hermes2d/examples/euler.rst
    src/hermes2d/examples/heat-transfer.rst
    src/hermes2d/examples/helmholtz.rst
    src/hermes2d/examples/linear-elasticity.rst
    src/hermes2d/examples/maxwell.rst
    src/hermes2d/examples/navier-stokes.rst
    src/hermes2d/examples/nernst-planck.rst
    src/hermes2d/examples/neutronics.rst
    src/hermes2d/examples/poisson.rst
    src/hermes2d/examples/richards.rst
    src/hermes2d/examples/schroedinger.rst
    src/hermes2d/examples/thermoelasticity.rst
    src/hermes2d/examples/wave-equation.rst
    src/hermes2d/examples/miscellaneous.rst
    
7. Hermes3D
-----------

.. toctree::
    :maxdepth: 1

    src/hermes3d/mesh.rst
    src/hermes3d/paraview.rst
    src/hermes3d/benchmarks.rst
    src/hermes3d/examples.rst

8. Hermes1D
-----------

.. toctree::
    :maxdepth: 1

    src/hermes1d/examples.rst



.. #####

    src/wrappers

    Indices and Tables
    ==================

    * :ref:`genindex`
    * :ref:`modindex`

    .. * :ref:`glossary`

    * :ref:`search`

    .. _Hermes: http://www.hpfem.org/hermes
    .. _FEMhub: http://www.hpfem.org/femhub
    .. _Agros2D: http://www.hpfem.org/agros2d
    .. _hp-FEM: http://www.hpfem.org
