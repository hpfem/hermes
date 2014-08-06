Mumps
-----

.. _MUMPS home page: `<http://graal.ens-lyon.fr/MUMPS/index.php>`_.

Linux
~~~~~

Using standard Debian packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For sequential version, install package `libmumps-seq-dev`. 
For parallel version, install `libmumps-dev` (this will require some changes in Hermes CMake configuration). In Ubuntu, you can use the Synaptic package manager for that, or type::

   sudo apt-get install libmumps-seq-dev

for the sequential version and
::

   sudo apt-get install libmumps-dev
   
for the parallel one.

Now go to the directory with Hermes. Create the file CMake.vars with the
following line (or append to the existing one)::

  set(WITH_MUMPS YES)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make

Find more about :ref:`ref-usage-mumps`.

Windows MSVC
~~~~~~~~~~~~

Installation of MUMPS using MSVC is rather easy:
  
  - preparation
  
    - download MUMPS from `<http://mumps.enseeiht.fr/MUMPS_4.10.0.tar.gz`_ (if the link does not work, look for 4.10 version of MUMPS)
    - download WinMUMPS utility from `<http://sourceforge.net/projects/winmumps/>`_.
    - download a Fortran compiler (e.g. `<http://software.intel.com/en-us/intel-fortran-studio-xe-evaluation-options)>`_.
    - download BLAS (Debug/Release, static/dynamic, 32-bit/64-bit as you like) from `<http://icl.cs.utk.edu/lapack-for-windows/lapack/index.html#libraries>`_.
    - you have to have Visual Studio version >= 2008
    - you have to have Python 2.6 or 2.7 available

  - installation
  
    - copy the downloaded BLAS library to 'bin' and 'lib' dependencies directory respectively. Please not that the name of the static library (static part of the dynamic library) should be either blas.lib, or libblas.lib
    - unzip MUMPS and WinMUMPS
    - execute "python winmumps_generator.py" from WinMUMPS providing the options needed (type -h for help)
    - this will generate several C++ project files (.vcxproj) and several Fortran project files (.vfproj)
    - open one of them and add the rest to the solution
    - build the solution in your desired configuration (Win32 or x64, Debug or Release)
    - this will place the libraries into MUMPS_4.10.0\lib\Debug\Win32\ (or alternatively for other configurations)
    - copy those to 'lib' dependencies directory
    - copy the contents of MUMPS_4.10.0\include\ to 'include' dependencies directory

Mac OS
~~~~~~

Help needed!

.. _ref-usage-mumps:

Using MUMPS in Hermes
~~~~~~~~~~~~~~~~~~~~~

After the installation has been completed, you may select  ``SOLVER_MUMPS`` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial` (in tutorial found on Hermes website), or use it just to solve a standalone matrix problem `Ax = b`.
