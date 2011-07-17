=============================
Welcome to the Hermes Project
=============================

This is a basic README file with generic information relevant for both Hermes Common and Hermes2D.


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


Installation
===========
::

    $ make install
		
The default installation paths are::

			/usr/local on Linux
			"C:\Program Files" on Win32
			"C:\Program Files (x86)" on Win64
