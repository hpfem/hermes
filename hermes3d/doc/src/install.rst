============
Installation
============

Linux Users
-----------

Compilation
~~~~~~~~~~~


Normally, you will need to satisfy some prerequisities before you can compile and install Hermes3D.
These are:

* cmake
* make
* C++ compiler
* Fortran compiler
* Judy
* BLAS
* UMFPACK

Clone the Git repository, configure and build:

.. code-block:: bash

   $ git clone http://hpfem.org/git/hermes3d.git 
   $ cmake . 
   $ make

If you have more than one CPU, you can use "make -jN" where N is the 
number of CPUs of your computer. 

Tests
~~~~~

To execute all test, do:

.. code-block:: bash

   $ make test

Note that some of the tests take a long time to finish. To execute the short 
running tests, do:

.. code-block:: bash

   $ make test-quick

Debugging with Eclipse
~~~~~~~~~~~~~~~~~~~~~~

To use eclipse as debugger, in the root folder of the project:

::

    mkdir eclipse_build
    cd eclipse_build
    cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../

In Eclipse:

    - Import project using Menu File->Import
    - Select General->Existing projects into workspace:
    - Browse where your build tree is and select the root build tree directory.
    - Keep "Copy projects into workspace" unchecked.

Install Hermes3D
~~~~~~~~~~~~~~~~

In order to install Hermes3D into different location, you need to change 
``CMAKE_INSTALL_PREFIX`` variable in CMakeCache.txt or issue command:

.. code-block:: bash

   $ cmake -DCMAKE_INSTALL_PREFIX=~/usr/local .
   $ make

To install the binary issue command:

.. code-block:: bash

   $ sudo make install		# if your install path is system-wide

Customization
~~~~~~~~~~~~~

Hermes3D is quite modular and can be build with several options. Configuration 
can be done via ``CMake.vars`` file that has to be placed in the build 
directory. Look into ``CMake.vars.example`` file for the file format.

- Element types

  * ``WITH_TETRA`` -- enable/dislable teterahedral elements.
  * ``WITH_HEX`` -- enable/dislable hexahedral elements.

- Library type:

  * ``REAL`` -- build real version of the library (scalar type is ``double``)
  * ``COMPLEX`` -- build complex version of the library (scalar type is ``std::complex<double>``)

- Modules

  * ``WITH_UMFPACK`` -- build with support for UMFPACK solver.
  * ``WITH_PETSC`` -- build with support for PETSc solver.
  * ``WITH_PARDISO`` -- build with support for PARDISO solver.
  * ``WITH_MUMPS`` -- build with support for MUMPS solver.
  * ``WITH_TRILINOS`` -- build with Trilinos support (AztecOO, Epetra, ML, Ifpack, NOX, Amesos).
  * ``WITH_METIS`` -- build with METIS support.
  * ``WITH_OPENMP`` -- build with OpenMP support.
  * ``WITH_MPI`` -- build with MPI support (currently no effect).

  If you have problems with CMake not finding your packages, you might want to check 
  ``cmake/FindXYZ.cmake`` files for further details and configuration options.

- Misc

  * ``DEBUG`` -- build debugging version of the library.
  * ``DEBUG_ORDER`` -- use the maximal integration order for integral evaluation.
  * ``WITH_TESTS`` -- build the tests to check that Hermes3D is doing what it is supposed to.
  * ``DEV_TESTS`` -- build developers tests. It is not recommended for normal users, these tests
    take very long time to finish (approx. weeks)
  * ``ADDITIONAL_LIBS`` -- the list of additional libraries that you need to link the binary files
    against in case you have some more requirements to fulfill. For example, if your PETSc is
    compiled with X11 support you need to link against X11 libs. Then you list all the libraries
    here.
  * ``OUTPUT_DIR`` -- set this to a directory that will contain the output files like solutions,
    convergence graphs, etc. Setting this to ``"."`` will write these files into the current
    directory.

For example, you can turn on and off various components to build, just edit 
the CMake.vars file and add the following:

::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

For development, edit the CMake.vars file and add the following:

:: 
    set(DEBUG YES)	# to compile debug versions
    set(RELEASE YES)	# to compile debug versions

Then issue command:

.. code-block:: bash

   $ make debug		# to build debug versions
   $ make release	# to build release versions

The CMake.vars can also be easily generated from a script (e.g. a debian/rules 
file) by issuing command:

.. code-block:: bash
    $ python -c 'print "set(H3D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

Another way how to configure the build process is to use ``ccmake``. Then the most of the options
listed above can be set via textual UI. The paths where your packages sit have to be specified in
``CMake.vars`` the same way as mentioned above.

Mac OS X Users
--------------

Compilation
~~~~~~~~~~~
**Step 1**: Make sure you have XCode installed. This should be on the 
installation disks which came with your Mac. XCode contains the GNU 
compilers, make and many other things which are required to build Hermes3D.

**Step 2**: Download and install MacPython version 2.6 using the disk image for
your version of OSX at http://www.python.org/download/releases/2.6.5/.
You will already have a version of Python which gets installed with
your operating system, but it will probably be out of date. Once this
is installed, go to the Python 2.6 directory which will be in your
Applications folder and double click the 'Update Shell
Profile.command' script to run it. This will update your system to use
the latest version of Python.

**Step 3**: Install the following libraries and applications: judy, Suitesparse,
glew, cmake, git. If you don't already have these on your Mac, then
the easiest way to get them is to use MacPorts (which is an
application which allows you to easily install and manage UNIX
libraries and applications on your Mac) by doing the following:

  (a) Download and install MacPorts from
      http://www.macports.org/install.php.
  (b) Do 'sudo port install judy suitesparse glew Lapack'.
  (c) If you don't already have git installed, do
      'sudo port install git'.
  (d) If you don't already have cmake installed, do
      'sudo port install cmake'.

**Step 4**: Get the Hermes3D source code. Change to the directory where you want
to download the Hermes3D source and clone the git repository by doing
'git clone http://hpfem.org/git/hermes3d.git'.

**Step 5**: Configure and build Hermes by doing 'cd hermes3d/ && cmake .
&& make'.
If you have more than one CPU, you can use 'make -jN' where N is the
number of CPUs of your computer. To set the location where Hermes2D
will be installed, pass the -DCMAKE_INSTALL_PREFIX=<your location>
flag to cmake (i.e. to install in /usr/local, replace the cmake
command above with 'cmake -DCMAKE_INSTALL_PREFIX=/usr/local .').

**Step 6**: To execute all tests, do 'make test'. Note that some of the tests 
can take a long time to finish. To just execute the short running tests,
do 'make test-quick'.

**Step 7**: Install Hermes3D by doing 'make install'.

Windows Cygwin Users
--------------------

Compilation
~~~~~~~~~~~

Download and install the Linux emulator Cygwin from `here <http://www.cygwin.com/>`_ (the small icon in the top-right corner). While running setup.exe, 
you need to install cmake, gcc4, gfortran, git, gitk, Lapack, libX11-devel, libXext-devel, libXt-devel, libXt, libXext, make, m4, openssl-devel, perl,
python, wget, xextproto.

Then download, unpack, and build Hermes3D as in Linux:

.. code-block:: bash

    $ git clone http://hpfem.org/git/hermes3d.git
    $ cd hermes3d
    $ cmake .
    $ make

For more details go to the Linux section above.

Notes
-----

* To build documentation, you will need to install the following packages:

   - Doxygen
   - breathe
   - sphinx (0.6.1 works)
   - TeX
   - dvipng

  ``make html`` builds the documentation in html format; there is also ``make latex`` which builds
  the TeX files with the documentation.

* When building the complex version of Hermes3D with PETSc support, you will need PETSc build with
  C++ support (i.e. ``--with-clanguage=C++`` when building PETSc)

* **Warning:** if you try to build the release version of Hermes3D on 64 bit machine using gcc 4.4.x
  you will get an error. This seems to be a problem in gcc, not in Hermes3D. To workaround this
  issue, build the debug version, if you cannot use different compiler or different version of gcc.   
 
* **Trilinos** can be build using our howto_.


.. _howto: http://hpfem.org/main/howto/howto-trilinos.html
