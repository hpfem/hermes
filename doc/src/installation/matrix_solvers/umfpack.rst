UMFpack
-------

.. _UMFPack home page: http://www.cise.ufl.edu/research/sparse/umfpack/
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/UMF-UserGuide.pdf

Linux
~~~~~

Using standard Debian packages
``````````````````````````````
Install the `libsuitesparse-metis-3.1.0` and `libsuitesparse-dev` packages.
In Ubuntu 9.10 (Karmic) or newer you can use the Synaptic package manager for that, or type::

    sudo apt-get install libsuitesparse-metis-3.1.0 libsuitesparse-dev

Now go to the directory with Hermes. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_UMFPACK YES)

and execute::

  rm CMakeCache.txt
  cmake .
  make
  
Find more about :ref:`ref-usage-umfpack`.

Using the special Hermes/Femhub package
```````````````````````````````````````

.. popup:: '#umfpack-matrix-solver'
   ../../../_static/clapper.png

Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory:

.. sourcecode::
   .

   wget https://github.com/downloads/hpfem/solvers/umfpack-5.5.1.spkg --no-check-certificate
   tar -jxvf umfpack-5.5.1.spkg
   rm umfpack-5.5.1.spkg
   cd umfpack-5.5.1

.. latexcode::
   .

   wget https://github.com/downloads/hpfem/solvers/umfpack-5.5.1.spkg 
   --no-check-certificate
   tar -jxvf umfpack-5.5.1.spkg
   rm umfpack-5.5.1.spkg
   cd umfpack-5.5.1

In order to install the library into say ``~/solvers/umfpack`` (you may choose any
path you like, provided that you have write access to it; the target directory 
will be created if it doesn't exist), type now into the terminal::

  ./standalone-install ~/solvers/umfpack

For advanced configuration possibilities, please read the `manual`_ or visit the 
`UMFPack home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space or 
just remove the object files by executing the following commands

::

  cd UFconfig; make clean
  cd AMD     ; make clean
  cd UMFPACK ; make clean

Now go to the directory with Hermes. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_UMFPACK YES)
  set(UMFPACK_ROOT ~/solvers/umfpack) #(or your own installation destination)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make

Find more about :ref:`ref-usage-umfpack`.

Windows (Cygwin, MinGW, MSVC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

http://matrixprogramming.com/2008/03/umfpack

Mac OS
~~~~~~

http://mywiki-science.wikispaces.com/UMFPACK

.. _ref-usage-umfpack:

Using UMFPACK in Hermes
~~~~~~~~~~~~~~~~~~~~~~~

After the installation has been completed, you may select ``SOLVER_UMFPACK`` as the matrix solver for your finite element problem,
as detailed in the `Poisson tutorial <http://http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`__, or use
it just to solve a standalone matrix problem :math:`Ax = b` as in the 
`Using Matrix Solvers tutorial <http://http://hpfem.org/hermes/doc/src/hermes2d/P08-miscellaneous/35-matrix-solvers.html>`__.
