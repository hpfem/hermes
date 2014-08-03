Trilinos
--------

.. _Trilinos home page: `<http://trilinos.sandia.gov/>`_.

Linux
~~~~~

Using standard Debian packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install packages `libtrilinos` and `libtrilinos-dev`. In Ubuntu 6.06 (Dapper)
or newer, you can use the Synaptic package manager for that, or type::

   sudo apt-get install libtrilinos libtrilinos-dev

Now go to the directory with Hermes. Create the file CMake.vars with the
following line (or append to the existing one)::

  set(WITH_TRILINOS YES)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make

Find more about :ref:`ref-usage-trilinos`.

Windows
~~~~~~~
First of all - to build Trilinos, one needs LAPACK (CLAPACK) (see the optional package in the library installation documentation).

| Download the sources for the latest version from the `Trilinos page <http://trilinos.sandia.gov>`__ and unpack them in some temporary directory.
| 
| Go to the Trilinos source directory.
| 
| In the following, replace {TPL_LAPACK_LIBRARIES}, {TPL_BLAS_LIBRARIES} with the full path to your lapack.lib and blas.lib without any quotes.
| Also, replace {CMAKE_INSTALL_PREFIX} with either your dependency root, or any other folder where you want to install Trilinos packages.

::

    mkdir build_dir
    cd build_dir
    cmake \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG \
     -D TPL_LAPACK_LIBRARIES:FILEPATH=d:\\hpfem\\hermes\\dependencies\\lib\\lapack.lib \
     -D TPL_BLAS_LIBRARIES:FILEPATH=d:\\hpfem\\hermes\\dependencies\\lib\\blas.lib \
     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
     -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
     -D Trilinos_ENABLE_Teuchos:BOOL=ON \
     -D Teuchos_ENABLE_TESTS:STRING=OFF \
     -D Teuchos_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_Epetra:BOOL=ON \
     -D Epetra_ENABLE_TESTS:STRING=OFF \
     -D Epetra_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
     -D EpetraExt_ENABLE_TESTS:STRING=OFF \
     -D EpetraExt_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_AztecOO:BOOL=ON \
     -D AztecOO_ENABLE_TESTS:STRING=OFF \
     -D AztecOO_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_Ifpack:BOOL=ON \
     -D Ifpack_ENABLE_TESTS:STRING=OFF \
     -D Ifpack_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_ML:BOOL=ON \
     -D ML_ENABLE_TESTS:STRING=OFF \
     -D ML_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_OpenMP:BOOL=ON \
     -D Trilinos_ENABLE_Amesos:BOOL=ON \
     -D Amesos_ENABLE_TESTS:STRING=OFF \
     -D Amesos_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_NOX:BOOL=ON \
     -D NOX_ENABLE_TESTS:STRING=OFF \
     -D NOX_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_Anasazi:BOOL=ON \
     -D Anasazi_ENABLE_TESTS:STRING=OFF \
     -D Anasazi_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_Komplex:BOOL=ON \
     -D Komplex_ENABLE_TESTS:STRING=OFF \
     -D Komplex_ENABLE_EXAMPLES:STRING=OFF \
     -D Trilinos_ENABLE_TESTS:BOOL=OFF \
     -D DART_TESTING_TIMEOUT:STRING=600 \
     -D CMAKE_INSTALL_PREFIX:STRING=/d/hpfem/hermes/dependencies \	 
     ..
	
| Build the Trilinos solution.
| Go up to the Trilinos source directory and edit the cmake_install.cmake file. Change::

	SET(CMAKE_INSTALL_CONFIG_NAME "Release")
	
for::

	SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
	
Install Trilinos into the path specified by the {CMAKE_INSTALL_PREFIX} variable by running::

	cmake -P cmake_install.cmake 
	
You may also need to create a dummy file "unistd.h" in the include folder under dependencies folder. This header is
not present in certain versions of Microsoft C library.
Go to the directory with Hermes. Add the following lines into CMake.vars::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT {CMAKE_INSTALL_PREFIX}) 
	
again, replace {CMAKE_INSTALL_PREFIX} with the folder where you installed Trilinos.

Find more about :ref:`ref-usage-trilinos`.

MAC OS
~~~~~~

In preparation.

.. _ref-usage-trilinos:

Using TRILINOS in Hermes
~~~~~~~~~~~~~~~~~~~~~~~~

You may now select either ``SOLVER_AMESOS`` as the direct matrix solver or 
``SOLVER_AZTECOO`` as the iterative matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P08-miscellaneous/35-matrix-solvers.html>`__.
Note that Trilinos is also required for using the advanced nonlinear solver ``NOX`` (see e.g. the 
`Trilinos - Nonlinear tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P07-trilinos/02-trilinos-nonlinear.html>`__).
