Trilinos
--------

.. _Trilinos home page: http://trilinos.sandia.gov/
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/Trilinos10.6Tutorial.pdf

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

Using the special Hermes/Femhub package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. popup:: '#trilinos-matrix-solver'
   ../../../_static/clapper.png

Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory:

.. sourcecode::
   .
  
   wget https://github.com/downloads/hpfem/solvers/trilinos-10.6.2.spkg --no-check-certificate
   tar -jxvf trilinos-10.6.2.spkg
   rm trilinos-10.6.2.spkg
   cd trilinos-10.6.2

.. latexcode::
   .
  
   wget https://github.com/downloads/hpfem/solvers/trilinos-10.6.2.spkg
   --no-check-certificate
   tar -jxvf trilinos-10.6.2.spkg
   rm trilinos-10.6.2.spkg
   cd trilinos-10.6.2

You can enable or disable various components of Trilinos in the script ``standalone_install``.
For example, if you do not want to enable PyTrilinos, change the corresponding line to::

   -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \

In order to install the library into say ``~/solvers/trilinos`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

  ./standalone-install ~/solvers/trilinos

For advanced configuration possibilities, please read the `manual`_ or visit the 
`Trilinos home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space or 
just remove the object files by executing the following command

::

  cd bin-dir; make clean

Now go to the directory with Hermes. Create the file CMake.vars with the following lines (or append to the existing one)::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT ~/solvers/trilinos) #(or your modified CMAKE_INSTALL_PREFIX variable)

Then execute::

    rm CMakeCache.txt
    cmake .
    make
    
Find more about :ref:`ref-usage-trilinos`.

Build from source
^^^^^^^^^^^^^^^^^

Download the sources for the latest version from the `Trilinos download page <http://trilinos.sandia.gov/download/trilinos-10.6.html>`__ and unpack them in some temporary directory. Go to the Trilinos source directory and issue the following commands there::

    mkdir build_dir
    cd build_dir
    cmake \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG \
     -D CMAKE_C_FLAGS:STRING="-fPIC -Wl,-V" \
     -D CMAKE_CXX_FLAGS:STRING="-fPIC -Wl,-V" \
     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
     -D BUILD_SHARED_LIBS:BOOL=ON \
     -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
     -D Trilinos_ENABLE_Teuchos:BOOL=ON \
     -D Trilinos_ENABLE_Epetra:BOOL=ON \
     -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
     -D Trilinos_ENABLE_AztecOO:BOOL=ON \
     -D Trilinos_ENABLE_Ifpack:BOOL=ON \
     -D Trilinos_ENABLE_ML:BOOL=ON \
     -D Trilinos_ENABLE_Amesos:BOOL=ON \
     -D Trilinos_ENABLE_NOX:BOOL=ON \
     -D Trilinos_ENABLE_Anasazi:BOOL=ON \
     -D Trilinos_ENABLE_Komplex:BOOL=ON \
     -D Trilinos_ENABLE_Claps:BOOL=ON \
     -D Trilinos_ENABLE_TESTS:BOOL=ON \
     -D Trilinos_ENABLE_MPI:BOOL=OFF \
     -D DART_TESTING_TIMEOUT:STRING=600 \
     -D CMAKE_INSTALL_PREFIX:STRING=~/solvers/trilinos \
     ..
  
    make
    sudo make install

(This installs the library into ~/solvers/trilinos directory. If you do not like this location, change the CMAKE_INSTALL_PREFIX variable to whatever you like.)

Now go to the directory with Hermes. Create the file CMake.vars with the following lines (or append to the existing one)::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT ~/solvers/trilinos) #(or your modified CMAKE_INSTALL_PREFIX variable)

Then execute::

    rm CMakeCache.txt
    cmake .
    make
    
Find more about :ref:`ref-usage-trilinos`.

Windows
~~~~~~~
| Download the sources for the latest version from the `Trilinos page <http://trilinos.sandia.gov/download/trilinos-10.6.html>`__ and unpack them in some temporary directory.
| 
| Go to the Trilinos source directory.
| 
| In the following, replace {CLAPACK_DIR} with the full path to your clapack-3.2.1-CMAKE directory (where you installed CLAPACK as a Hermes's dependency) without any quotes.
| Also, replace {CMAKE_INSTALL_PREFIX} with either your dependency root, or any other folder where you want to install Trilinos packages.::


    mkdir build_dir
    cd build_dir
    cmake \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG \
     -D CLAPACK_DIR:STRING={CLAPACK_DIR} \
     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
     -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
     -D Trilinos_ENABLE_Teuchos:BOOL=ON \
     -D Trilinos_ENABLE_Epetra:BOOL=ON \
     -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
     -D Trilinos_ENABLE_AztecOO:BOOL=ON \
     -D Trilinos_ENABLE_Ifpack:BOOL=ON \
     -D Trilinos_ENABLE_ML:BOOL=ON \
     -D Trilinos_ENABLE_Amesos:BOOL=ON \
     -D Trilinos_ENABLE_NOX:BOOL=ON \
     -D Trilinos_ENABLE_Anasazi:BOOL=ON \
     -D Trilinos_ENABLE_Komplex:BOOL=ON \
     -D Trilinos_ENABLE_Claps:BOOL=ON \
     -D Trilinos_ENABLE_TESTS:BOOL=ON \
     -D Trilinos_ENABLE_MPI:BOOL=OFF \
     -D DART_TESTING_TIMEOUT:STRING=600 \
     -D CMAKE_INSTALL_PREFIX:STRING={CMAKE_INSTALL_PREFIX} \	 
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
