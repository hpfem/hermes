Mumps
-----

.. _MUMPS home page: http://graal.ens-lyon.fr/MUMPS/index.php
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/MUMPS_4.9.2.pdf

Linux
~~~~~

Using standard Debian packages
``````````````````````````````
Install the libmumps-seq-dev and libmumps-seq-4.9.2 packages.
In Ubuntu you can use the Synaptic package manager for that, or type::

    sudo apt-get install libmumps-seq-dev libmumps-seq-4.9.2

Note that these packages are available since Ubuntu 10.10 (Maverick). If you have an older version of Ubuntu or want an up-to date version of the library, please follow the instructions below.

Using the special Hermes/Femhub package
```````````````````````````````````````
Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory::
  
  wget https://github.com/hpfem/solvers/raw/master/packages/mumps-4.9.2.spkg --no-check-certificate
  tar -jxvf mumps-4.9.2.spkg
  rm mumps-4.9.2.spkg
  cd mumps-4.9.2

In order to install the library into say ``~/solvers/mumps`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

  SPKG_LOCAL=~/solvers/mumps ./spkg-install

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

You may now select ``SOLVER_MUMPS`` as the matrix solver for your finite element problem, as detailed
in the `Poisson tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-5/matrix_solvers.html>`__.

Windows MSVC
~~~~~~~~~~~~

http://matrixprogramming.com/2010/05/mumps

Mac OS
~~~~~~

Help needed!
