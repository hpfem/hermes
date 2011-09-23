====================
Hermes Documentation
====================

User Documentation in HTML
--------------------------

Before building User Documentation, install the Python Sphinx package::

    sudo apt-get install python-sphinx

The user documentation can be found in the directory doc/. Type "make html" there 
to build it. The HTML pages can then be displayed by typing::

    firefox _build/html

or using another web browser. 

User Documentation in PDF
-------------------------

In the directory doc/ type "make latex" and you will be instructed how to build 
the PDF::

    Build finished; the LaTeX files are in _build/latex.
    Run `make all-pdf' or `make all-ps' in that directory to run these through (pdf)latex.



Developer Documentation
-----------------------

In order to view developer documentation, install Doxygen::

    sudo apt-get install doxygen

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
about Hermes and adaptive *hp*-FEM can be found in the `publications section <http://hpfem.org/publications/>`_
at hpfem.org. The most recent list is probably the one 
on `Pavel Solin's publications page <http://hpfem.org/~pavel/public/papers.html>`_.


