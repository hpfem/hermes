============
About Hermes
============

Hermes is a free C++/Python library for rapid development of
adaptive *hp*-FEM and *hp*-DG solvers for partial differential equations (PDE)
and multiphysics PDE systems. The developer team includes the 
`hp-FEM group <http://hpfem.org/>`_ at the `University of Nevada, Reno <http://www.unr.edu>`_ 
and their `collaborators <http://git.hpfem.org/git/hermes.git/tree/HEAD:/AUTHORS>`_ 
from numerous places around the globe.

The standard way to use Hermes is to write short C++ user programs, but for 
those who prefer to use a graphical interface, there is an 
`interactive GUI Agros2D <http://hpfem.org/agros2d/>`_. We also provide 
an `interactive online lab <http://nb.femhub.org/>`_ where
the user can compute with Hermes and other FEM codes in `FEMhub <http://femhub.org>`_ 
via any web browser without having to install anything (CPU time is on us). 

Although Hermes is much younger than major FEM packages, it is loaded with 
unique technology and its user base is growing fast. We hope that you will 
enjoy the software and that you will find this documentation useful. Let us know if 
you find mistakes or with any improvement suggestions. Anyone who contributes
a patch becomes automatically a 
`co-author <http://git.hpfem.org/git/hermes.git/tree/HEAD:/AUTHORS>`_ of the code.

The library is available under the GPL license (Version 2, 1991). 

User and Developer Documentation
--------------------------------

User documentation (tutorial, benchmarks, examples) can be found in
the directory 'doc/'. Type 'make html' there to build it. The HTML pages are then
in _build/html.

In the following, you can replace "2d" with "3d" or "1d" to obtain links
for Hermes3D and Hermes1D, respectively.

To compile the C++ reference manual for Hermes2D, go to 'hermes2d/doc.cpp/'. There
type 'doxygen hermes2d.lib-real.doxyfile' to build references for the 
real version, or 'doxygen hermes2d.lib-cplx.doxyfile' to build refs for the 
complex version. The HTML files are then in 'h2d-real/html/index.html' and
'h2d-cplx/html/index.html', respectively. The Doxygen documentation for the 
real and complex version is also 
available online at http://hpfem.org/hermes2d/doc.cpp/h2d-real/html/index.html
and http://hpfem.org/hermes2d/doc.cpp/h2d-cplx/html/index.html, respectively.

.. toctree::
   :maxdepth: 2

    math-background
    web-access
    citing-hermes


