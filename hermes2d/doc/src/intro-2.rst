============
Installation
============

Linux Users
-----------

Compilation
~~~~~~~~~~~

If you are using a Debian-based system, install the required libraries first:

:: 

    apt-get install cmake g++ gfortran libjudy-dev freeglut3-dev libsuitesparse-dev libglew1.5-dev python-dev python-numpy python-scipy cython python-matplotlib

(Note: cmake has to be at least version 2.6 or later, matplotlib has to be at
least 0.98.5.2 or higher.)

Clone the Git repository, configure and build:

::
  
    git clone http://git.hpfem.org/git/hermes2d.git
    cd hermes2d/
    cmake .
    make

If you have more than one CPU, you can use "make -jN" where N is
the number of CPUs of your computer.

Tests
~~~~~

To execute all tests, do:

::

    make test

Note that some of the tests take a long time to finish. To just execute the
short running tests, do:

::

    make test-quick

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by:

::

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

If you are on OS X, you have to disable GLUT as the glut library is not easily
installable on OS X. To do so, just put the following line into your
CMake.vars:

::

    set(WITH_GLUT NO)


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


Install Hermes2D
~~~~~~~~~~~~~~~~

::

    cmake -DCMAKE_INSTALL_PREFIX=~/usr .
    make
    make install

Mac OS X Users
--------------

Compilation
~~~~~~~~~~~

**Step 1**: Make sure you have XCode installed. This should be on the installation 
disks which came with your Mac. XCode contains the GNU compilers, make 
and many other things which are required to build Hermes2D.

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
  (b) Do 'sudo port install judy suitesparse glew'.
  (c) If you don't already have git installed, do 
      'sudo port install git'.
  (d) If you don't already have cmake installed, do 
      'sudo port install cmake'.

**Step 4**: Get the Hermes2D source code. Change to the directory where you want 
to download the Hermes2D source and clone the git repository by doing 
'git clone http://git.hpfem.org/git/hermes2d.git'.

**Step 5**: Configure and build Hermes by doing 'cd hermes2d/ && cmake . 
&& make'.
If you have more than one CPU, you can use 'make -jN' where N is the 
number of CPUs of your computer. To set the location where Hermes2D 
will be installed, pass the -DCMAKE_INSTALL_PREFIX=<your location> 
flag to cmake (i.e. to install in /usr/local, replace the cmake 
command above with 'cmake -DCMAKE_INSTALL_PREFIX=/usr/local .').

**Step 6**: To execute all tests, do 'make test'. Note that some of the tests can 
take a long time to finish. To just execute the short running tests, 
do 'make test-quick'.

**Step 7**: Install Hermes2D by doing 'make install'.

Tests
~~~~~

To execute all tests, do:

::
 
    make test

Note that some of the tests take a long time to finish. To just execute the
short running tests, do:

::

    make test-quick


More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by:

::

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars


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

Compilation
~~~~~~~~~~~

Download and install the Linux emulator Cygwin from `here <http://www.cygwin.com/>`_ (the small icon in the top-right corner). While running setup.exe, you need to install 

cmake, gcc4, gfortran, git, gitk, libX11-devel, libXext-devel, libXt-devel, libXt, libXext, make, m4, openssl-devel, perl, 
python, wget, xextproto.

Then download, unpack, and build FEMhub as in Linux:

::

    git clone http://git.hpfem.org/git/hermes2d.git
    cd hermes2d
    cmake .
    make

For more details go to the Linux section above.

Windows MSVC Users
------------------

This section describes how to build and use Hermes2D in Microsoft Visual C++ 2008 (Express Edition). 
These instructions should probably work even for older versions of MS Visual C++ up to version 2003.

Known limitations and issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - Stand-alone viewers do not compile.
 - Trilinos not supported.
 - Python not supported.

Building Hermes2D
~~~~~~~~~~~~~~~~~

 In order to build the library and examples, you need to:

 - Prepare dependecy libraries, see 'Dependency Check-list'.
 - Copy a file 'my_hermes_root\\MSVC2008\\CMake.vars' to 'my_hermes_root'. The file contains settings for the projekct.
 - Modify the file 'my_hermes_root\\CMake.vars'. Only the first line has to be modified, the rest of lines depends on project settings. Please, follow comments in the file. For example, you 
   could set the first line as:  

   ::

       set(DEP_ROOT "../dependencies")

 - In the directory 'my_hermes_root', to create project files by running CMAKE from a command prompt:

   ::

       cmake . -G "Visual Studio 9 2008"  # MSVC2008 user
       cmake . -G "Visual Studio 10"      # MSVC2010 user

   If you have Cygwin installed, your might have an error "Coulld not create named generator Visual Studio 10". This is because your 
   cmake path is contaminated by Cygwin's cmake. Try to use absoluted path for windows cmake.exe. 

 - Open a SLN file 'my_hermes_root\\hermes2d.sln' and build Hermes2D. Actually, this step is not necessary if you plan to use Hermes2D in your projects rather than explore tutorials.

Configuration options
~~~~~~~~~~~~~~~~~~~~~

 Hermes2D is configured through preprocessor directives. Directives are generated by CMAKE and your settings might be overriden by CMAKE. The directives are:

  - NOGLUT : excludes GLUT-dependant parts. This replaces viewers with an empty implementation that does nothing if invoked. If used, a library 'freeglut.lib' does not need to be linked. 

  - ENABLE_VIEWER_GUI : enables GUI for viewers. Currently, only 'ScalarView' support GUI because this is an experimental feature. This directive is mutually exclusive with NOGLUT. If used, a library 'AntTweakBar.lib' does not need to be linked.

Using Hermes2D
~~~~~~~~~~~~~~
 
In order to used Hermes2D in your project, you need to do following steps. Steps has 5, 6, and 7 to be repeated for every configuration, i.e., Debug, Release. Except the step 7b, this can be done easily by setting the drop-down Configuration to 'All configurations' in the Project Property dialog.

  - Prepare Hermes2D to be buildable by MSVC, see 'Building Hermes2D'.
  - Create your project in MSVC. Set the project to be empty Win32 console project.
  - Add either 'my_hermes2d_root\\src\\hermes2d-real.vcproj' or 'my_hermes2d_root\\src\\hermes2d-cplx.vcproj' project to your solution (<right click on solution>\\Add\\Existing Project...).
  - Set that your project depends on hermes2d-real/cplx project (<right click on your project>\\Project Dependences...).
  - Add directories 'my_hermes2d_directory\\src' and 'dependencies\\include' to additional include directories (<right click on your project>\\Properties\\Configuration Properties\\C/C++\\Additional Include Directories).
  - Add directories 'dependencies\\lib' to additional library directories (<right click on your project>\\Properties\\Configuration Properties\\Linker\\Additional Library Directories).
  - Deny warnings that are not indicating anything dangerous:

    - Avoid warnings about STL in DLL by denying a warning 4251 (<right click on your project>\\Properties\\Configuration Properties\\C/C++\\Advanced\\Disable Specific Warnings, enter 4251).
    - Avoid warnings about standard functions that are not safe (<right click on your project>\\Properties\\Configuration Properties\\C/C++\\Preprocessor\\Preprocessor Definitions, add _CRT_SECURE_NO_WARNINGS).
 
Dependency check-list
~~~~~~~~~~~~~~~~~~~~~

This list works for 32-bit version of Hermes2D. If you intend to cross-compile 64-bit version, you have to cross-compile all libraries. As the first step, create a  directory structure.
	
  - In order to create the structure, execute the following command: 

    ::

        'prepare_dep_dir.bat'. 

    Be sure to include a directory 'dependecies\\bin' into 'PATH' environment variable.
  - All Hermes2D project files assumes that dependency libraries are available in a fixed directory structure. The root of this structure has to have the same parent as does Hermes2D director, i.e., if 'C:\\my_work\\hermes2d\\' is a root of the Hermes2D directory, then 'C:\\my_work\\dependecies\\' is a root of the dependency directory. Subdirectories are:

    - dependencies\\include: Header files (\*.h) of dependency libraries.
    - dependencies\\lib: Library files (\*.lib) of dependency libraries.   
    - dependencies\\bin: Binary modules (\*.dll) of dependency libraries. Be sure to include a directory 'dependecies\\bin' into 'PATH' environment variable.


  - Python(ver-2.6.5)

    - Download Python 2.6.5 windows installer (http://www.python.org/ftp/python/2.6.5/python-2.6.5.msi) and run it. 
    - If you Python installation path is 'C:\\Python26', add it to System Path variable. You may need to reboot to have the Path variable updated. 
    - Open a command prompt, execute python, and you will see something like:

      ::

          python.exe
          Python 2.6.5 (r265:79096, Mar 19 2010, 21:48:26) [MSC v.1500 32 bit (Intel)] on win32
          >>>

      If you have Cygwin installed again, make sure you are running the Windows Python. 
 
    - Numpy(ver-1.4.1)

      - Download Numpy-1.4.1-py2.6-superpack-python2.6.exe (http://sourceforge.net/projects/numpy/files/) and run it.
      - Open a command prompt, excute following command from python shell:

        ::

            >>> import numpy
            >>> print numpy.get_include()
            C:\Pytthon26\lib\site-packages\numpy\core\include
            >>>

    - Cython(ver-0.12.1)

      - Download Cython (http://www.cython.org/#download), and unpack it. 
      - Open a command prompt, got to 'my_Cython_root', and excute following: 

        ::

            python.exe setup.py install
	
  - JUDY(ver-1.0.5)

    - Download judy (http://sourceforge.net/projects/judy/) and upack it.
    - Open a command prompt with MSVC variables set up: either use a command prompt in MSVC start menu or start a command prompt and execute VCVARS.BAT from the MSVC directory.
    - Switch to a directory containing JUDY sources, e.g., 'my_judy_root\\src'.
    - Compile JUDY with 'build.bat': this creates Judy.dll and Judy.lib.
    - Copy 'Judy.dll', 'Judy.h', and 'Judy.lib' to 'bin', 'include', and 'lib' dependecy directories respectively.

  - PTHREAD(2.8.0)

    - Download pthread binaries version 2.8.0 (ftp://sourceware.org/pub/pthreads-win32/prebuilt-dll-2-8-0-release/).
    - Copy 'lib\\pthreadVCE2.dll', 'include\\\*.h' and 'lib\\pthreadVCE2.lib' to 'bin', 'include', and 'lib' dependecy directories respectively.

  - UMFPACK(ver-5.4.0)

    - UFConfig(ver-3.4.0):

      - Download UFconfig source file package (http://www.cise.ufl.edu/research/sparse/UFconfig/), and unpack it. 
      - Copy UFconfig.h to 'include' dependecy directory.

    - AMD(ver-2.2.0):

      - Download AMD source file package (http://www.cise.ufl.edu/research/sparse/amd/).
      - Unpack source file into a directory that has the same parent as a directory where you unpacked UFconfig.
      - Copy the file 'my_hermes2d_root\\MSVC2008\\AMD.nmake' to a directory 'my_amd_directory\\Lib'.
      - Run MSVC command prompt similar as in the case of Judy and switch to 'my_amd_directory\\Lib'.
      - Compile AMD using 'nmake -f AMD.nmake'.
      - Copy 'Include\\amd.h', 'Include\\amd_internal.h', and 'Lib\\libamd.lib' to 'include', and 'lib' dependecy directories respectively.

    - Download UMFPACK source file package (http://www.cise.ufl.edu/research/sparse/umfpack/).
    - Unpack source file into a directory that has the same parent as a directory where you unpacked UFconfig.
    - Copy the file 'my_hermes2d_root\\MSVC2008\\UMFPACK.nmake' to the 'my_umfpack_root\\Lib' directory.
    - Run MSVC command prompt similar as in the case of Judy and switch to 'my_umfpack_root\\Lib'.
    - Compile UMFPACK using 'nmake -f UMFPACK.nmake'. Linking might take some time, please, be patient.
    - Copy 'libumfpack.dll', all include files, and 'libumfpack.lib' to 'bin', 'include', and 'lib' dependecy directories, respectively.

  - CMAKE

    - Download CMAKE (http://www.cmake.org/cmake/resources/software.html) version 2.6.4 source, and cmake 2.8.1 binary. 
    - Since 2.8.1 came out after MSVC10, you need to download (http://www.cmake.org/files/v2.8/CMakeVS10FindMake.cmake) and 
      replace your_cmake_2.8.1_root\\share\\cmake-2.8\\Modules\\CMakeVS10FindMake.cmake, If you are using MSVC10. 
    - Double click cmake-gui (2.8.1), choose your_cmake_2.6.4_root as source directory, and your_cmake_2.6.4_tmp as build directory. 
      click configure button, wait for a while.  After configuration is done, change CMAKE_INSTALL_PREFIX as: C:\\my_work\\dependencies 
      so that cmake is installed under the 'dependencies\\bin' and accessible from every location. 
    - Now, go to your_cmake_2.6.4_tmp (build dir), you will find a CMAKE.sln file asociated with MSVC. Open it with MSVC08/10, and 
      find project "INSTALL", right click, choose "Build". 
 
  - OpenGL support (optional)

    - If a directive NOGLUT is used, this step and all its substeps can be skipped.
    - FREEGLUT 

      - Download freeglut 2.4.0 (http://freeglut.sourceforge.net/) and unpack it.
      - Open the your_freeglut_2.4.0_root\\freeglut.DSP file in MSVC08, MSVC08 will convert file into a newer format, 
          i.e., SLN and VCPROJ (MSVC10 user could open the converted file freeglut.sln).
      - Compile Debug or Release version. Debug version is recommended in a case of debugging.
      - Copy 'freeglut.dll', 'freeglut.h', and 'freeglut.lib' to 'bin', 'include\\GL', and 'lib' dependency directories, respectively/.
  
  - GLEW

    - Download glew 1.5.4 (http://glew.sourceforge.net/) and unpack it.
    - Open a DSP file 'my_glew_root\\builds\\vc6\\gle_shred.dsp' and let MSVC08 to 
        convert it (MSVC10 user could open the converted file glew_shared.sln).
    - Switch to 'Release' version, and  build a project 'glew_shared': this will create DLL file.
    - Copy 'my_glew_root\\bin\\glew32.dll', 'my_glew_root\\include\\GL\\\*.h', and 'my_glew_root\\lib\\glew32.lib' 
        to 'bin', 'include\\GL', and 'lib' dependency directories respectively.
 	
  - AntTweakBar (optional)

    - If a directive ENABLE_VIEWER_GUI is *not* used, this step can be skipped.
    - Download a modified version 1.1.3 of AntTweakView (http://hpfem.org/downloads/AntTweakBar.1.1.3.modified.tar.gz) and unpack it. 
    - Open SLN file in MSVC and compile it.
    - Copy 'AntTweakBar.dll', 'AntTweakBar.h', and 'AntTweakBar.lib' to 'bin', 'include', and 'lib' dependency directories respectively.
	
  - ExodusII (optional)

    - If a directive WITH_EXODUSII is *not* used, this step including all sub-steps can be skipped.
	
    - Zlib

      - Download sources of version 1.2.3 (http://sourceforge.net/projects/libpng/files/) and unpack them.
      - Open 'my_zlib_root/projects/visualc6/zlib.dsw' (Visual C++ 6 Solution File) in MSVC08 and let MSVC to convert it and save the .sln file 
          (MSVC10 user can open the .sln file).
      - Switch a configuration to 'Release DLL' in Configuration Manager. 
      - Build project 'zlib': this will create DLL/LIB files in 'my_zlib_root/projects/visual6/Win32_DLL_Release'.
      - Copy 'zlib1.dll', 'zlib.h/zconf.h', and 'zlib1.lib' to 'bin', 'include', and 'lib' dependency directories respectively.
 
    - HDF5

      - Download sources of version 1.8.x (ftp://ftp.hdfgroup.org/HDF5/hdf5-1.8.0/src/) and unpack them. 
      - Since SLIB is not used, comment out a line '#define H5_HAVE_FILTER_SZIP 1' in the header file 'my_hdf5_root\\windows\\src\\H5pubconf.h'
      - Copy the file 'my_hdf5_root\\windows\\src\\H5pubconf.h' to the directory 'my_hdf5_root\\src\\'
      - Run MSVC Command Prompt and switch to a directory 'my_hdf5_root\\windows\\proj'
      - Set variable HDF5_EXT_ZLIB to 'my_dependencies\\lib\\zlib1.lib', by issusing the following:

        ::

            set HDF5_EXT_ZLIB="C:\my_hermes2d_root\dependencies\lib\zlib1.lib


      - If SLIB is used, set variable HDF5_EXT_SLIB similarly as:

        ::

            set HDF5_EXT_SLIB="C:\my_hermes2d_root\dependencies\lib\slib.lib

      - To open SLN file in MSVC by issusing the following in the command prompot, and let MSVC to convert files: 

        ::

            VCExpress.exe all\all.sln

      - Switch a configuration to 'Release'
      - Build project 'hdf5_hldll': this will create DLL/LIB files in 'my_hdf5_root\\proj\\hdf5_hldll\\Release\\' and 'my_hdf5_root\\proj\\hdf5dll\\Release\\'
      - Copy 'hdf5dll.dll' and 'hdf5dll.lib' to 'bin' and 'lib' dependency directories respectively
      - Copy 'hdf5_hldll.dll' and 'hdf5_hldll.lib' to 'bin' and 'lib' dependency directories respectively
      - Currently, only MSVC08 is supported under Vista. But MSVC08/10 should be supported under Windows XP. 

    - NetCDF

      - Download sources of version 4.0.1 (http://www.unidata.ucar.edu/downloads/netcdf/netcdf-4_0_1/index.jsp) and unpack them.
      - Open a SLN file 'my_netcfd_root\\win32\\NET\\netcdf.sln'.
      - Switch to 'Release' version.
      - In properties of the project 'netcdf'. 

        - Add paths 'my_hdf5_root\\src\\' and 'my_hdf5_root\\hl\\src' to 'C/C++ -> Additional Include Directories'
        - Add a path 'dependencies\\lib\\' to 'Linker -> Additional Library Directories'

      - Build project 'netcdf': this will create DLL/LIB files in 'my_netcdf_root/win32/NET/Release'
      - Copy 'netcdf.dll' and 'netcdf.lib' to 'bin' and 'lib' dependency directories respectively
      - Copy 'my_netcdf_root\\libsrc4\\netcdf.h' to 'include' dependency directory

    - ExodusII

      - Download sources of version 4.9.3 (http://sourceforge.net/projects/exodusii/) and unpack 'exodusii'
      - Add the following line to the file 'my_exodusii_root\\CMakeLists.txt' as:

        ::

            PROJECT(Exodusii)
            SET(NETCDF_INCLUDE_DIR "my_netcdf_root/libsrc4")    # add this line; 

        be sure to use a slash '/' instead of a backslash '\\'. 

      - Generate MSVC project files using CMAKE in command prompt as:

        ::

            cmake . -G "Visual Studio 9 2008"    # MSVC2008 user 
            cmake . -G "Visual Studio 10"        # MSVC2010 user 

        If you have Cygwin installed, make sure that you are using the windows version of cmake. 

      - Open a SLN file 'my_exodusii_root/ExodusII.sln' in MSVC08/10
      - Switch to 'Release' version
      - Build a project 'exoIIv2c': this will create a LIB file in 'my_exodusii_root\\cbind\\Release'
      - Copy 'exoIIv2c.lib' to 'lib' dependency directory structure
      - Copy 'my_exodusii_root\\cbind\\include\\exodusII.h and exodusII_ext.h' to 'include' dependency directory
	

=============================
How to Submit the First Patch
=============================

The following is an embarrassingly trivial git primer
whose objective is to show you how to create and send 
your first patch without losing much time and good humor. 
We begin with cloning the Hermes2D git repository and 
continue through setting 
up the .gitconfig file, creating a new branch, committing 
changes, and generating patches. A good reference for 
further reading is given at the end. 

Clone the Hermes2D Git Repository
---------------------------------

To clone the repository, type

::

    git clone http://git.hpfem.org/git/hermes2d.git

This will create a new directory hermes2d/ with a copy 
of the entire Hermes2D git repository. Before doing anything 
else, you may want to build Hermes2D to make sure that 
everything works:

::

    cd hermes2d/
    cmake .
    make

The list of prerequisites and installation instructions 
for various platforms can be found 
`here <http://hpfem.org/hermes2d/doc/src/intro-2.html>`_.

Create the .gitconfig File
--------------------------

The .gitconfig file can be used to define your identity
for git as well as useful abbreviations. In hermes2d/
type 

::

    cd .. 

Then adjust and save the following as "~/.gitconfig":

::

    [user]
	    name = Pavel Solin
	    email = solin.pavel@gmail.com

    [core]
	    editor = vim

    [color]
	    ui = true
    [color "branch"]
	    current = yellow reverse
	    local = yellow
	    remote = green
    [color "diff"]
	    meta = yellow bold
	    frag = magenta bold
	    old = red bold
	    new = green bold
	    whitespace = red reverse
    [color "status"]
	    added = yellow
	    changed = green
	    untracked = cyan
    [core]
	    whitespace=fix,-indent-with-non-tab,trailing-space,cr-at-eol

    [alias]
	    st = status
	    ci = commit
	    br = branch
	    co = checkout
	    df = diff
	    lg = log -p

Create a Local Branch
---------------------

Type 

::

    cd hermes2d/ 

You can get an overview of existing branches by typing 

::

    git branch 

This will show you something like this:

  .. image:: img/intro/terminal-git.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

If this is your first time, then you will see
just the master branch with the star next to it,
which tells you that there are no other branches.

If you want to make any changes to the source files, then 
it is a good idea to always create a new branch for it. 
This is done by typing

::

    git co -b test-1

where test-1 is the name of your new local branch. Now you 
can do any changes you like and you do not have to be afraid
of damaging your master branch. HOWEVER, you always must 
commit your changes as described below. 
Unless you commit your changes, git does not 
know that they belong to your local branch. Then you are not 
able to switch branches, and YOU ARE IN TROUBLE!

Make and Commit Your Changes
----------------------------

The best way to get your first patch in 
is to look into the source files in the 
directory src/, find any function that 
is missing a comment, and fix this. 
Say that this file was src/file.cpp.
After adding the comment there, type:

::

    git add src/file.cpp
    git commit

The latter command will invoke a basic text editor 
where you will be asked to enter a one-line comment
describing your changes. If you decide to skip this 
and commit with an empty line, your commit will not 
be accepted. 

Create and Send a Patch
-----------------------

You are almost there! Just type 

::

    git format-patch -1

and a new text file starting with three zeros will be 
created. This is a "patch". The parameter '-1' in there
means that you want only the last commit included in 
the patch. If you typed '-2', git would include the last 
two commits, etc. 

Last, send an email with the patch to the mailing 
list hermes2d@googlegroups.com, begin the subject 
line with saying "[PATCH] ...", and attach the 
text file with the patch to your email. Someone
will be with you shortly!

Change to Master and Update the Repository
------------------------------------------

Before changing to a different branch, type 

::

    git st

This stands for 'git status'. You will see 
something like this:

  .. image:: img/intro/terminal-git-2.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

The green font tells you that git has the latest 
version of the file. All modified files in red 
need to be added using "git add". It is a good
idea to go through the untracked files too, in case
that you wish to add some of them as well. 
Related to the sample screenshot above, after 
typing 

::

    git add src/intro-2.rst
    git st

you will see

  .. image:: img/intro/terminal-git-3.png
   :align: center
   :width: 600
   :alt: Terminal screenshot

Now you can proceed with "git commit" as described above. 
After the commit, you can switch to the master branch:

::

    git co master

This brings you to the point where you can 
return to the beginning of this short
tutorial, and start working on a new change.

To update your master to the latest state of
the repository, just type:

::

    git remote add origin http://git.hpfem.org/git/hermes2d.git

This tells git where to download the git repository from
(needs to be done just the first time). Then type

::

    git pull origin master

Special Note on Sphinx Docs
---------------------------

The Sphinx documentation you are just reading is also 
part of the Hermes2D git repository and can be developed
in the same way as source files of Hermes2D. This very 
file can be found in doc/src/intro-2.rst. After making 
any changes to Sphinx docs, type 

::

    make html

in the doc/ directory. This will launch 
a build process that will take a few seconds. 
After it finishes, type

::

    firefox _build/html

This will open a new tab in your Firefox where you will
see something like 

  .. image:: img/intro/firefox.png
   :align: center
   :width: 600
   :alt: Firefox screenshot

Click on the link "index.html" and you should see
the local form of your Sphinx docs that include your 
changes. 

Further Reading and Video
-------------------------

Git is very powerful and we covered just a tiny part of 
it. After the above works for you, please
read more about git in `Pro Git <http://progit.org/book/>`_.

Also watch this `YouTube video <http://www.youtube.com/watch?v=OFkgSjRnay4>`_
by Scott Chacon.

GitHub
------

You should also learn how to upload
your local branch to `GitHub <http://github.com/>`_
instead of sending a patch, since this makes the
work with your changes easier. 

Good luck and let us know if you think 
that this document could be improved!


 







