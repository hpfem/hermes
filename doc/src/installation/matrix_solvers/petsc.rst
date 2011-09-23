PETSc
-----

.. _PETSc home page: http://www.mcs.anl.gov/petsc/
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/petsc.pdf

Linux
~~~~~

Using standard Debian packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install packages `libpetsc3.1` and `libpetsc3.1-dev`. 
Newer version may be available. In Ubuntu 6.06 (Dapper)
or newer, you can use the Synaptic package manager for that, or type::

   sudo apt-get install libpetsc3.1 libpetsc3.1-dev

Now go to the directory with Hermes. Create the file CMake.vars with the
following line (or append to the existing one)::

  set(WITH_PETSC YES)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make

Find more about :ref:`ref-usage-petsc`.

Using the special Hermes/Femhub package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. popup:: '#petsc-matrix-solver'
   ../../../_static/clapper.png

.. only:: latex

    `Tutorial Video <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/videos.html#petsc-matrix-solver>`_. 

Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory:

.. sourcecode::
   .

   wget https://github.com/downloads/hpfem/solvers/petsc-3.1-p7.spkg --no-check-certificate
   tar -jxvf petsc-3.1-p7.spkg
   rm petsc-3.1-p7.spkg
   cd petsc-3.1-p7

.. latexcode::
   .

   wget https://github.com/downloads/hpfem/solvers/petsc-3.1-p7.spkg
   --no-check-certificate
   tar -jxvf petsc-3.1-p7.spkg
   rm petsc-3.1-p7.spkg
   cd petsc-3.1-p7

In order to install the library into say ``~/solvers/petsc`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

  ./standalone-install ~/solvers/petsc

For advanced configuration possibilities, please read the `manual`_ or visit the 
`PETSc home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space.

Now go to the directory with Hermes. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_PETSC YES)
  set(PETSC_ROOT ~/solvers/petsc) #(or your own installation destination)
  set(PETSC_ARCH linux-cxx)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make
  
Find more about :ref:`ref-usage-petsc`.

Windows MSVC
~~~~~~~~~~~~

http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html#Windows

Mac OS
~~~~~~

http://petsc.darwinports.com/

.. _ref-usage-petsc:

Using PETSC in Hermes
~~~~~~~~~~~~~~~~~~~~~

You may now select ``SOLVER_PETSC`` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P08-miscellaneous/35-matrix-solvers.html>`__.
