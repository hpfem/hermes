Mesh
~~~~
First part one needs to handle is the computational mesh, typically the following would be used::

    // Shared pointers are used for easier memory handling.
    MeshSharedPtr mesh(new Mesh);
    // Native Hermes mesh format.
    MeshReaderH2D mloader;
    mloader.load("domain.mesh", mesh);
    
    // XML Hermes mesh format.
    MeshReaderH2DXML mloader;  
    mloader.load("domain.xml", mesh);
    
More about meshes can be found in the 'hermes-tutorial' documentation, section 'A-linear', chapter '01-mesh'.