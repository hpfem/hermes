=============================
Welcome to the Hermes Project
=============================

This is a basic README file with generic information relevant for 
all Hermes1D, Herme2D and Hermes3D. The directories hermes1d/,
hermes2d/ and hermes3d/ contain their own README.rst files with 
specific details. 


Copyright
=========

Copyright (c) 2010 hp-FEM group at the University of Nevada,
Reno (UNR). Email: hpfem@unr.edu, home page: http://hpfem.org/.


User Documentation
==================

User Documentation is available online at http://hpfem.org/hermes/doc/index.html.


Download and Installation
=========================

Installation instructions for Linux, Mac OS X, Windows Cygwin,
and Windows MSVC are part of User Documentation. 


Generating and Viewing User Documentation Offline
=================================================

Linux
-----

Source files of the (Sphinx) user documentation are in the
directory doc/. In order to compile the user documentation,
you need to install Sphinx. Then follow these steps::

    cd doc
    make html
    firefox _build/html

Windows + MSVC
--------------

This is a sequence of steps which will install Sphinx and which
will allow to generate user documentation. The steps assume that you
have already installed Microsoft Visual Studio 9 (MSVC) Express
Edition (or any higher edition) and you have a copy of Hermes sources.

All commands, which are marked with a keyword 'prompt:', are executed
in a command prompt opened in the step 3. Search for all mentioned
application through Google since, usually, the first link is the right
one.

#. Download and install python 2.6
#. Add paths 'my_python_path\' and 'my_python_path\Scripts' to
   the enviromental variable PATH.
#. Open a command prompt with MSVC variables:
   Search for 'Visual Studio 2008 Command Prompt' in the start menu.
#. Download and install setuptools 0.6c11
#. Install Sphinx using setuptools (Internet access required)
   prompt: easy_install -U Sphinx
#. Go to a folder doc folder of Hermes source tree
   prompt: cd my_hermes\doc
#. Run NMAKE requesting HTML version
   prompt: nmake html
#. View the documentation using a file
   'my_hermes2d\doc\_build\index.html'


Generating and Viewing Developer Documentation Offline
======================================================

Linux
-----

Source files of the (Doxygen) developer documentation are
kept in the directories hermes1d/doc.cpp, hermes2d/doc.cpp
and hermes3d/doc.cpp. In order to build them you need
to install Doxygen. If you are in hermesXd/,
do::

    cd doc.cpp/
    doxygen hermesXd.lib-real.doxyfile
    doxygen hermesXd.lib-cplx.doxyfile

This will generate documentation for the real and complex
version of the library, respectively. To view the docs,
type::

    firefox hXd-real/html/index.html
    firefox hXd-cplx/html/index.html


Compilation
===========

::

    $ cmake .
    $ make

If you have cmake text-based UI installed, you can do::

    $ ccmake .
      press C
      customize your build
      press C
      press G
    $ make

Hermes build configuration scripts will look for required libraries on default
include and system paths. If you packages are installed elsewhere, you need to
specify their paths. Consult either online documetation (see the link above) or
cmake/FindXYX.cmake files for exact names. These values have to be put into
CMake.var file located in your build directory.

Intel C Compiler
----------------

To use Intel C compiler::

    $ export CC=/path/to/icc
    $ export CXX=/path/to/icpc
    $ cmake .
    $ make

NOTE: version 10.0.026 did not work for us (some compatibility issues with
STL), 10.1.022 works ok
