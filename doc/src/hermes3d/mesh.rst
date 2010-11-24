.. _examples-doc:

Mesh Data Format
================

Hermes3D can read meshes in a native format and in 
`ExodusII format <http://sourceforge.net/projects/exodusii/>`_.
For the former, use

::

    // Load the mesh file.
    Mesh mesh;
    H3DReader mloader;
    mloader.load("domain.mesh", &mesh);

at the beginning of your main.cpp file. 
To load an ExodusII mesh file, use ``ExodusIIReader`` class instead
of "H3DReader". 

Below is a sample mesh file containing a single hex element. Note the 
ordering of vertices in elements and on faces. Also note 
that indices start from 1 and that elements do not 
have material markers as in Hermes2D. We are going to 
resolve both these incompatibilities soon::

    # vertices
    8
    -1 -1  -1
     1 -1  -1
     1  1  -1
    -1  1  -1

    -1 -1   1
     1 -1   1 
     1  1   1
    -1  1   1

    # tetras
    0

    # hexes
    1
    1 2 3 4 5 6 7 8

    # prisms
    0 

    # tris
    0 

    # quads
    6
    1 2 3 4		5
    1 2 6 5		3
    2 3 7 6		2
    3 4 8 7		4
    4 1 5 8		1
    5 6 7 8		6

In the quads section, the last number on each line 
is a boundary marker as the reader expects.






