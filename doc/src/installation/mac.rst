Mac OS
======

Download and compilation
~~~~~~~~~~~~~~~~~~~~~~~~

**Known issues**: Hermes has built-in OpenGL visualization based on FreeGlut, but this package 
comes with certain installation difficulties. If you encounter Glut-related problems, set 
H2D_WITH_GLUT to NO in Cmake.vars, build Hermes without Glut, and use VTK output for visualization. 

**Step 1**: Make sure you have XCode installed. This should be on the installation 
disks which came with your Mac. XCode contains the GNU compilers, make 
and many other things which are required to build Hermes.

**Step 3**: Install the following libraries and applications: Suitesparse, 
glew, cmake, git. If you don't already have these on your Mac, then 
the easiest way to get them is to use MacPorts (which is an 
application which allows you to easily install and manage UNIX 
libraries and applications on your Mac) by doing the following:

  (a) Download and install MacPorts from 
      http://www.macports.org/install.php.
  (b) Do 'sudo port install suitesparse glew'.
  (c) If you don't already have git installed, do 
      'sudo port install git'.
  (d) If you don't already have cmake installed, do 
      'sudo port install cmake'.

**Step 4**: Get the Hermes source code as described at the beginning of the Linux section
above. Change to the directory where you want 
to download the Hermes source and clone the git repository either
from the hpfem.org server::

    git clone http://git.hpfem.org/git/hermes.git

or from Github::

    git clone git://github.com/hpfem/hermes.git

These two repositories are synchronized. For more advanced users we recommend 
to create a free account at Github (if you do not have one yet), fork the 
Hermes repository, and then clone your Github copy of Hermes to your local computer. 
This will establish links between your local copy and the master repository, and 
you’ll become part of the Hermes network at Github.

**Step 5**: Configure and build Hermes by changing dir to 'hermes/', 
and then typing 'cmake .' and 'make'.
If you have more than one CPU, you can use 'make -jN' where N is the 
number of CPUs of your computer. To set the location where Hermes 
will be installed, pass the -DCMAKE_INSTALL_PREFIX=<your location> 
flag to cmake (i.e. to install in /usr/local, replace the cmake 
command above with 'cmake -DCMAKE_INSTALL_PREFIX=/usr/local .').

**Step 6**: Install Hermes by doing 'make install'.

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following::

    set(WITH_EXAMPLES NO)
    
(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by:

.. sourcecode::
    .

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

.. latexcode::
    .

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)
    \nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars


For development, it is good to say (in global CMake.vars)::

    set(DEBUG YES) to compile debug versions
    set(RELEASE YES) to compile release versions

Then type::

    make debug    (to build debug versions)
    make release  (to build release versions)
