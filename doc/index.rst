.. Hermes2D documentation master file, created by
   sphinx-quickstart on Fri Oct 30 19:23:04 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================
Hermes Documentation
====================

Thank you for your interest in Hermes! 

This document is organized as follows: 
Section 1 contains general information about Hermes, its math background,
and interactive web accessibility. You will also learn how to use Git and 
Github, and how to contribute to the project if interested. Section 2 will 
walk you through its installation on various platforms and related issues. 
The tutorial itself begins with Section 3 that explains in many incremental 
steps how Hermes2D is used. Sections 4 and 5 contain examples and benchmarks 
for Hermes3D and Hermes1D (their usage is practically identical to Hermes2D).

This document is under continuous development and certainly it is not perfect. 
If you find bugs, typos, dead links or such, help us improve it by reporting them
through one of the mailing lists for 
`Hermes1D <http://groups.google.com/group/hermes1d/>`_,
`Hermes2D <http://groups.google.com/group/hermes2d/>`_, or
`Hermes3D <http://groups.google.com/group/hermes3d/>`_. 

Thank you in advance for helping us improve Hermes!

1. Introduction
---------------

.. toctree::
    :maxdepth: 1

    src/about-hermes
    src/math-background
    src/web-access
    src/citing-hermes
    src/first_pull_request

2. Installation
---------------

.. toctree::
    :maxdepth: 1

    src/installation/linux
    src/installation/mac
    src/installation/win-cygwin
    src/installation/win-msvc
    src/installation/matrix_solvers

3. Hermes2D
-----------

.. toctree::
    :maxdepth: 1

    src/hermes2d/linear
    src/hermes2d/nonlinear
    src/hermes2d/timedep
    src/hermes2d/linear-adapt
    src/hermes2d/nonlinear-adapt
    src/hermes2d/timedep-adapt
    src/hermes2d/eigen
    src/hermes2d/fvm-and-dg
    src/hermes2d/trilinos
    src/hermes2d/miscellaneous
    src/hermes2d/benchmarks
    src/hermes2d/nist
    src/hermes2d/examples
    
4. Hermes3D
-----------

.. toctree::
    :maxdepth: 1

    src/hermes3d/mesh.rst
    src/hermes3d/benchmarks.rst
    src/hermes3d/examples.rst

5. Hermes1D
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
