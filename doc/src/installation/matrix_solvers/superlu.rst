SuperLU
--------

Hermes currently supports two versions of the SuperLU library - the sequential
one and the multithreaded one. Support for the MPI version will be added in the 
future. Please visit `<http://crd.lbl.gov/~xiaoye/SuperLU/>`_ for more information about the
library.

Linux
~~~~~

Using standard Debian packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the `libsuperlu3` and `libsuperlu3-dev` packages. In Ubuntu 6.06 (Dapper)
or newer, you can use the Synaptic package manager for that, or type::

  sudo apt-get install libsuperlu3 libsuperlu3-dev 
  
Now go to the directory with Hermes. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_SUPERLU YES)
  set(SUPERLU_ROOT ~/solvers/superlu_mt) #(or your own installation destination)
  set(SUPERLU_MT   NO)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make
  
Find more about :ref:`ref-usage-superlu`.

Windows MSVC
~~~~~~~~~~~~

http://crd.lbl.gov/~xiaoye/SuperLU/faq.html

MAC OS
~~~~~~

http://www.bleedingmind.com/index.php/2010/07/31/compiling-superlu-on-os-x/

.. _ref-usage-superlu:

Using SUPERLU in Hermes
~~~~~~~~~~~~~~~~~~~~~~~

You may now select ``SOLVER_SUPERLU`` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://http://hpfem.org/hermes/doc/src/hermes2d/P08-miscellaneous/35-matrix-solvers.html>`__.
