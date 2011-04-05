=============================================
Installation of ExodusII and NetCDF libraries
=============================================

ExodusII and NetCDF libraries are necessary to read Cubit files in Hermes.


Installing ExodusII
-------------------

To install the ExodusII library, insert the following commands in your command-line terminal 
(be root user for system wide installation):

Step 1)

:: 

    cd /

Step 2)

64bit:

::

    wget http://hpfem.org/downloads/exodusii-bin-amd64-4.81.tar.gz

32bit:

::

    wget http://hpfem.org/downloads/exodusii-bin-i686-4.81.tar.gz

Step 3)

::

    tar xvzf exodusii-bin-*.tar.gz

Step 4)

::

    rm exodusii-bin-*.tar.gz

The latter steps will unpack the binary package of ExodusII into /opt/packages/exodusii.  
If you do not like the location, you may change it to whatever you like, just remember to also adjust the "CMake" variables accordingly (see below).

To use ExodusII support in Hermes2d or 3d: - clone the Hermes repo (if you haven't done so yet), and 
in your "CMake.vars" file (found in your Hermes directory) add the following lines (AFTER you are done installing BOTH [ExodusII and NetCDF] packages)::

    SET(WITH_EXODUSII YES)
    SET(EXODUSII_ROOT /opt/packages/exodusii)
    SET(NETCDF_ROOT /opt/packages/netcdf)

After you are done adding these lines to your "CMake.vars" file, go back to your command-line 
terminal and type::

    cmake .
    make

Congratulations (assuming everything worked out) you just installed the ExodusII library!


Installing NetCDF
-----------------

To install the NetCDF library, insert the following commands in your command-line terminal
(be root user for system wide installation):

Step 1)

::

    cd /

Step 2)

64bit:

::

    wget http://hpfem.org/downloads/netcdf-bin-amd64-4.0.1.tar.gz

32bit:

::

    wget http://hpfem.org/downloads/netcdf-bin-i686-4.0.1.tar.gz

Step 3)

::

    tar xvzf netcdf-bin-*.tar.gz

Step 4)

::

    rm netcdf-bin-*.tar.gz

The latter steps will unpack the binary package of NetCDF into /opt/packages/netcdf.
If you do not like the location, you may change it to whatever you like, just remember to also adjust the "CMake" variables accordingly (see below).

To use ExodusII support in Hermes2d or 3d: - clone the Hermes repo (if you haven't done so yet), and
in your "CMake.vars" file (found in your Hermes directory) add the following lines (AFTER you are done installing BOTH [ExodusII and NetCDF] packages)::

    SET(WITH_EXODUSII YES)
    SET(EXODUSII_ROOT /opt/packages/exodusii)
    SET(NETCDF_ROOT /opt/packages/netcdf)

After you are done adding these lines to your "CMake.vars" file, go back to your command-line
terminal and type::

    cmake .
    make

Congratulations (assuming everything worked out) you just installed the NetCDF library!

Note: The binary packages were prepared on Ubuntu 9.10 (64bit and 32bit); they may or may not work on other
64bit and 32bit machines. 




































 







