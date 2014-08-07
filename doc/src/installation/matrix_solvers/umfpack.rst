UMFpack
-------

.. _UMFPack home page: `<http://www.cise.ufl.edu/research/sparse/umfpack/>`_.

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

and execute

[OPTIONAL] If you want to use fast saving / loading of Hermes entities, install

  - BSON
  
    - Clone the BSON Mongo driver git repository from https://github.com/l-korous/mongo-c-driver.git (if you don't know how, here is a tip: `Getting a Git Repository <http://git-scm.com/book/en/Git-Basics-Getting-a-Git-Repository>`_)
    - Compile and install using 'make install'

[OPTIONAL] For thread caching memory allocator from Google, see
    
  - TCMalloc
    
      - Get TCMalloc from the SVN repository at `<http://code.google.com/p/gperftools/source/checkout>`_.
      - Make & install::

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

After the installation has been completed, you may select ``SOLVER_UMFPACK`` as the matrix solver for your finite element problems::

  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_UMFPACK);
