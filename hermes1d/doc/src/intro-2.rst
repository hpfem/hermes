============
Installation
============

Linux Users
-----------

Documentation
~~~~~~~~~~~~~

Hermes1D Sphinx documentation can be found in
the directory 'doc/'. Type 'make html' there to build it. The documentation is
available online at http://hpfem.org/hermes1d/doc/index.html.

Developer documentation can be compiled by running 'doxygen' in 'src/'.


Compilation
~~~~~~~~~~~

If you are using a Debian-based system, install the required libraries first:

:: 

    apt-get install cmake g++ python-dev python-numpy python-scipy cython python-matplotlib

(Note: cmake has to be at least version 2.6 or later, matplotlib has to be at
least 0.98.5.2 or higher.)

Clone the Git repository, configure and build:

::
  
    git clone http://hpfem.org/git/hermes1d.git
    cd hermes1d/
    cmake .
    make

If you have more than one CPU, you can use "make -jN" where N is
the number of CPUs of your computer.

Tests
~~~~~

To execute all tests, do:

::

    make test

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_PYTHON YES)
    set(WITH_DEBUG NO)
    set(WITH_RELEASE YES)
    set(WITH_TESTS NO)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

For development, it is good to say (in global CMake.vars):

::

    set(DEBUG YES) to compile debug versions
    set(RELEASE YES) to compile release versions

Then type:

::
 
    make debug    (to build debug versions)
    make release  (to build release versions)

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


Install Hermes1D
~~~~~~~~~~~~~~~~

::

    cmake -DCMAKE_INSTALL_PREFIX=~/usr .
    make
    make install

Mac OS X Users
--------------

Documentation
~~~~~~~~~~~~~

Hermes1D Sphinx documentation can be found in
the directory 'doc/'. Type 'make html' there to build it. The documentation is
available online at http://hpfem.org/hermes1d/doc/index.html.

Developer documentation can be compiled by running 'doxygen' in 'src/'.

Compilation
~~~~~~~~~~~

**Step 1**: Make sure you have XCode installed. This should be on the installation 
disks which came with your Mac. XCode contains the GNU compilers, make 
and many other things which are required to build Hermes1D.

**Step 2**: Download and install MacPython version 2.6 using the disk image for 
your version of OSX at http://www.python.org/download/releases/2.6.5/. 
You will already have a version of Python which gets installed with 
your operating system, but it will probably be out of date. Once this 
is installed, go to the Python 2.6 directory which will be in your 
Applications folder and double click the 'Update Shell 
Profile.command' script to run it. This will update your system to use 
the latest version of Python.

**Step 3**: Install the following libraries and applications: 
cmake, git. If you don't already have these on your Mac, then 
the easiest way to get them is to use MacPorts (which is an 
application which allows you to easily install and manage UNIX 
libraries and applications on your Mac) by doing the following:

  (a) Download and install MacPorts from 
      http://www.macports.org/install.php.
  (b) If you don't already have git installed, do 
      'sudo port install git'.
  (c) If you don't already have cmake installed, do 
      'sudo port install cmake'.

**Step 4**: Get the Hermes1D source code. Change to the directory where you want 
to download the Hermes1D source and clone the git repository by doing 
'git clone http://hpfem.org/git/hermes1d.git'.

**Step 5**: Configure and build Hermes by doing 'cd hermes1d/ && cmake . 
&& make'.
If you have more than one CPU, you can use 'make -jN' where N is the 
number of CPUs of your computer. To set the location where Hermes1D 
will be installed, pass the -DCMAKE_INSTALL_PREFIX=<your location> 
flag to cmake (i.e. to install in /usr/local, replace the cmake 
command above with 'cmake -DCMAKE_INSTALL_PREFIX=/usr/local .').

**Step 6**: To execute all tests, do 'make test'.

**Step 7**: Install Hermes1D by doing 'make install'.

Tests
~~~~~

To execute all tests, do:

::
 
    make test

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_PYTHON YES)
    set(WITH_DEBUG NO)
    set(WITH_RELEASE YES)
    set(WITH_TESTS NO)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

For development, it is good to say (in global CMake.vars):

::

    set(DEBUG YES) to compile debug versions
    set(RELEASE YES) to compile release versions

Then type:

::
 
    make debug    (to build debug versions)
    make release  (to build release versions)

Windows Cygwin Users
--------------------

Download and install the Linux emulator Cygwin from `here <http://www.cygwin.com/>`_ (the small icon in the top-right corner). While running setup.exe, you need to install 

cmake, gcc4, gfortran, git, gitk, make, m4, openssl-devel, perl, 
python, wget, xextproto.

Then download, unpack, and build Hermes1D as in Linux:

::

    git clone http://hpfem.org/git/hermes1d.git
    cd hermes2d
    cmake .
    make

For more details go to the Linux section above.

 







