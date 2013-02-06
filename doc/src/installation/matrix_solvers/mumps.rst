Mumps
-----

.. _MUMPS home page: http://graal.ens-lyon.fr/MUMPS/index.php
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/MUMPS_4.9.2.pdf

Linux
~~~~~

Using standard Debian packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For sequential version, install packages `libmumps-seq-4.9.2` and `libmumps-seq-dev`. 
For parallel version, install `libmumps-4.9.2` and `libmumps-dev`. Newer versions 
may be available. In Ubuntu 6.06 (Dapper)
or newer, you can use the Synaptic package manager for that, or type::

   sudo apt-get install libmumps-seq-4.9.2 libmumps-seq-dev

for the sequential version and
::

   sudo apt-get install libmumps-4.9.2 libmumps-dev
   
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

http://matrixprogramming.com/2010/05/mumps

Mac OS
~~~~~~

Help needed!

.. _ref-usage-mumps:

Using MUMPS in Hermes
~~~~~~~~~~~~~~~~~~~~~

After the installation has been completed, you may select  ``SOLVER_MUMPS`` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/P08-miscellaneous/35-matrix-solvers.html>`__.
