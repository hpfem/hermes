Mumps
-----

.. _MUMPS home page: http://graal.ens-lyon.fr/MUMPS/index.php
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/MUMPS_4.9.2.pdf

Linux
~~~~~

Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory::
  
  wget https://github.com/downloads/hpfem/solvers/mumps-4.9.2.spkg --no-check-certificate
  tar -jxvf mumps-4.9.2.spkg
  rm mumps-4.9.2.spkg
  cd mumps-4.9.2

In order to install the library into say ``~/solvers/mumps`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

  ./standalone-install ~/solvers/mumps

For advanced configuration possibilities, please read the `manual`_ or visit the
`MUMPS home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space or 
just remove the object files by running

::

  cd src
  make clean 

Now go to the directory with Hermes. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_MUMPS YES)
  set(MUMPS_ROOT ~/solvers/mumps) #(or your own installation destination)

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
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-5/matrix_solvers.html>`__.
