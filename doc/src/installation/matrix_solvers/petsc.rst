PETSc
-----
.. _PETSc home page: www.mcs.anl.gov/petsc/
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/petsc.pdf

Linux
~~~~~

Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory::
  
  wget https://github.com/hpfem/solvers/raw/master/packages/petsc-3.1-p5.spkg --no-check-certificate
  tar -jxvf petsc-3.1-p5.spkg
  rm petsc-3.1-p5.spkg
  cd petsc-3.1-p5

In order to install the library into say ``~/solvers/petsc`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

SPKG_LOCAL=~/solvers/petsc ./spkg-install

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

You may now select `SOLVER_PETSC` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-5/matrix_solvers.html>`__.

Windows MSVC
~~~~~~~~~~~~

http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html#Windows

Mac OS
~~~~~~

http://petsc.darwinports.com/

