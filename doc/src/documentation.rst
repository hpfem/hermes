====================
Hermes Documentation overview
====================

Building user documentation (this one) in HTML
----------------------------------------------

Before building User Documentation, install the Python Sphinx package::

    Linux: sudo apt-get install python-sphinx
    Windows: download Python (http://python.org), download setup tools (http://pypi.python.org/pypi/setuptools), adjust PATH env. variable

The user documentation can be found in the directory doc/. Type "make html" there 
to build it. The HTML pages can then be displayed by typing::

    firefox _build/html

or using another web browser. 


Building user documentation (this one) in PDF
---------------------------------------------

In the directory doc/ type "make latex" and you will be instructed how to build 
the PDF::

    Build finished; the LaTeX files are in _build/latex.
    Run `make all-pdf' or `make all-ps' in that directory to run these through (pdf)latex.


Developer Documentation (in Doxygen)
------------------------------------
The documentation is accessible online.

     `Hermes - common code <http://hpfem.org/~hermes/hermes/hermes_common/doc/html/index.html>`_

     `Hermes - 2D specific code <http://hpfem.org/~hermes/hermes/hermes2d/doc/html/index.html>`_


In order to build developers documentation, install Doxygen::

    Linux: sudo apt-get install doxygen
    Windows: http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.2-setup.exe

There are separate Doxygen files for hermes-common (dimension-independent functionality
such as matrix solvers) and for Hermes2D. To build the former, go to the directory 
hermes-common/ and type "doxygen Doxyfile". The HTML docs will be located in a new
subfolder doc/html/ and you can view them in any web browser, such as::

    firefox doc/html/index.html

To build Doxygen files for Hermes2D, go to the directory hermes2d/ and type again
"doxygen Doxyfile". The HTML docs will be located in the file doc/html/index.html.


Other Resources
---------------

See the section "Citing Hermes" in this document for a representative selection of 
books and scientific papers about Hermes. A more complete overview of publications 
about Hermes and adaptive *hp*-FEM can be found in the `publications section <http://hpfem.org/people/>`_
at hpfem.org. The most recent list is probably the one 
on `Pavel Solin's publications page <http://hpfem.org/~pavel/public/papers.html>`_.


