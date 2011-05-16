============
About Hermes
============

Hermes is a free C++/Python library for rapid development of
adaptive *hp*-FEM and *hp*-DG solvers for partial differential equations (PDE)
and multiphysics PDE systems. The development team includes the 
`hp-FEM group <http://hpfem.org/>`_ at the `University of Nevada, Reno <http://www.unr.edu>`_ 
and their `collaborators <http://git.hpfem.org/hermes.git/blob/HEAD:/AUTHORS>`_ 
from numerous places around the globe.

A standard way to use Hermes is to write short C++ user programs 
that use the functionality provided by the library, but for 
those who prefer to use a graphical interface, there is 
`Agros2D <http://hpfem.org/agros2d/>`_. 

.. image:: img/agros.png
   :scale: 50 %   
   :align: center 	
   :alt: Agros snapshot.

Hermes is loaded with modern finite element technology and its user base is 
growing fast. We hope that you will enjoy the software and that you will find 
this documentation useful. In any case please let us know if you find mistakes 
or if you can suggest improvements to this documentation or to Hermes itself.
Anyone who contributes with at least one patch becomes automatically a 
`co-author <http://git.hpfem.org/hermes.git/blob/HEAD:/AUTHORS>`_ of the code.

The library is available under the GPL license (Version 2, 1991). 

User and Developer Documentation
--------------------------------

This user documentation can be found in the directory doc/. Type "make html" there 
to build it. The HTML pages can then be displayed by typing "firefox _build/html",
"chromium-browser _build/html" or similar.

To compile the C++ reference manual for Hermes2D, go to 'hermes2d/doc.cpp/'. There
type "doxygen hermes2d.lib-real.doxyfile" to build references for the 
real version, or "doxygen hermes2d.lib-cplx.doxyfile" to build refs for the 
complex version. The HTML files are then in "h2d-real/html/index.html" and
"h2d-cplx/html/index.html", respectively. The Doxygen documentation for the 
real and complex version is also 
available online at http://hpfem.org/hermes2d/doc.cpp/h2d-real/html/index.html
and http://hpfem.org/hermes2d/doc.cpp/h2d-cplx/html/index.html, respectively.

.. toctree::
   :maxdepth: 2

    math-background
    citing-hermes


