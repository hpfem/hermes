Linux
-----

Download and compilation
~~~~~~~~~~~~~~~~~~~~~~~~

[NEW] You can download a package directly from **`Hermes launchpad repository<https://launchpad.net/~lukas-korous/+archive/ubuntu/hermes>`_**

The rest of the instructions here are for building Hermes from source.

If you are using a Debian-based system, install the (required) libraries first:

.. sourcecode::
    .

    apt-get install git git-core cmake g++ freeglut3-dev libsuitesparse-dev libglew-dev libxerces-c-dev xsdcxx libmatio-dev

If you want to use fast saving / loading of Hermes entities, install

  - BSON
  
    - Clone the BSON Mongo driver git repository from git@github.com:l-korous/mongo-c-driver.git (if you don't know how, here is a tip: `Getting a Git Repository <http://git-scm.com/book/en/Git-Basics-Getting-a-Git-Repository>`_)
    - Compile and install using 'make install'

For thread caching memory allocator from Google, see
    
  - TCMalloc
    
      - Get TCMalloc from the SVN repository at http://code.google.com/p/gperftools/source/checkout
      - Make & install
  
To obtain the source code, clone the Git repository from Github::
  
    git clone git://github.com/hpfem/hermes.git

These two repositories are synchronized. For more advanced users we recommend to 
create a free account at `Github <http://github.com>`_ (if you do not have one yet),
fork the `Hermes repository <http://github.com/hpfem/hermes>`_, and then clone your 
Github copy of Hermes to your local computer. This will establish links between
your local copy and the master repository, and you'll become part of the Hermes 
network at Github.

Once you have a local copy of the Hermes repository on your computer, change dir 
to hermes/. There you will find a CMake.vars.example.Linux file that looks like this::

    # LINUX
      # On linux, there should be no need to set up *_ROOT directories, in the default settings, they all point to /usr/local, as should be true on Debian systems.
      # We mainly support gcc and CLang compilers with C++11 support.
      
      # BASIC CONFIGURATION
      
      # Global
    # Generate static libs (instead of dynamic)
      set(HERMES_STATIC_LIBS NO)
      # Target path
      set(CMAKE_INSTALL_PREFIX "/usr/local")
      
      # Paths for compulsory dependencies
      set(XERCES_ROOT "/usr/local")
      set(XSD_ROOT "/usr/local")
      
      # HermesCommon
        
        # Release and debug versions
        set(HERMES_COMMON_DEBUG     YES)
        set(HERMES_COMMON_RELEASE   YES)
        ...


Copy this file to "CMake.vars" and set the variables according to your needs.
After that, type::

    cmake .
    make

If you have more than one CPU, you can use "make -jN" where N is
the number of CPUs of your computer.

Debugging with Eclipse
~~~~~~~~~~~~~~~~~~~~~~

To use eclipse as debugger, in the root folder of the project::

    mkdir eclipse_build
    cd eclipse_build
    cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../

In Eclipse:

    - Import project using Menu File->Import
    - Select General->Existing projects into workspace:
    - Browse where your build tree is and select the root build tree directory. 
    - Keep "Copy projects into workspace" unchecked.


Install Hermes
~~~~~~~~~~~~~~

::

    make install
