Trilinos
--------

Linux
~~~~~

Download the sources for the latest version from the `Trilinos page <http://trilinos.sandia.gov/download/trilinos-10.6.html>`__ and unpack them in some temporary directory. Go to the Trilinos source directory and issue the following commands there::

    mkdir build_dir
    cd build_dir
    cmake \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG \
     -D CMAKE_C_FLAGS:STRING="-fPIC -Wl,-V" \
     -D CMAKE_CXX_FLAGS:STRING="-fPIC -Wl,-V" \
     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
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

Go to the directory with hermes{1|2|3}d. Create the file CMake.vars with the following lines (or append to the existing one)::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT ~/solvers/trilinos) #(or your modified CMAKE_INSTALL_PREFIX variable)

Then execute::

    rm CMakeCache.txt
    cmake .
    make

Windows
~~~~~~~
| First, you need to install CLAPACK/CBLAS:
| 
| Download the file clapack-3.2.1-CMAKE.tgz from http://www.netlib.org/clapack/.
| 
| Use cmake to configure and build the debug version of clapack.
| 
| Copy '\\clapack-3.2.1-CMAKE\\BLAS\\SRC\\Debug\\blas.lib', '\\clapack-3.2.1-CMAKE\\F2CLIBS\\libf2c\\Debug\\libf2c.lib', and '\\clapack-3.2.1-CMAKE\\SRC\\Debug\\lapack.lib' to 'lib' dependency directory.
| 
| Copy the contains of '\\clapack-3.2.1-CMAKE\\INCLUDE\\' to 'include' dependency directory.
 
| Download the sources for the latest version from the `Trilinos page <http://trilinos.sandia.gov/download/trilinos-10.6.html>`__ and unpack them in some temporary directory.
| 
| Go to the Trilinos source directory.
| 
| In the following, replace {CLAPACK_DIR} with the full path to your clapack-3.2.1-CMAKE directory without any quotes.
Also, replace {CMAKE_INSTALL_PREFIX} with either your dependency root, or any other folder where you want to install Trilinos packages.::


    mkdir build_dir
    cd build_dir
    cmake \
     -D CMAKE_BUILD_TYPE:STRING=DEBUG \
     -D CLAPACK_DIR:STRING={CLAPACK_DIR} \
     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
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
Go up to the Trilinos source directory and edit the cmake_install.cmake file. Change::

	SET(CMAKE_INSTALL_CONFIG_NAME "Release")
	
for::

	SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
	
Install Trilinos into the path specified by the {CMAKE_INSTALL_PREFIX} variable by running::

	cmake -P cmake_install.cmake 
	
Go to the directory with hermes{1|2|3}d. Add the following lines into CMake.vars::

    set(WITH_TRILINOS YES)
    set(TRILINOS_ROOT {CMAKE_INSTALL_PREFIX}) 
	
again, replace {CMAKE_INSTALL_PREFIX} with the folder where you installed Trilinos.

MAC OS
~~~~~~

In preparation.
