=======================
Installation of Solvers
=======================

Mumps
-----

In preparation.

Pardiso
-------

In preparation.

PETSc
-----

In preparation.

Trilinos
--------

Linux
~~~~~

Download the sources from the `Trilinos page <http://trilinos.sandia.gov/download/trilinos-10.4.html>`__ and unpack them in some temporary directory. Go to the trilinos source directory::

    mkdir build_dir
    cd build_dir
    wget http://hpfem.org/main/howto/trilinos/do-configure
    chmod +x do-configure
    ./do-configure
    make
    sudo make install

(This places the installation into /opt/packages/trilinos dir. If you do not like this location, change the CMAKE_INSTALL_PREFIX variable to whatever you like.)

Go to the directory with hermes{1|2|3}d. Create the file CMake.vars::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT /opt/packages/trilinos) #(or your modified CMAKE_INSTALL_PREFIX variable)

And then execute::

    cmake .
    make

Windows
~~~~~~~

In preparation.

MAC OS
~~~~~~

In preparation.

UMFpack
-------

In preparation.






