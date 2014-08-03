First steps
-----------

After a successful compilation of Hermes library, you can install it by typing 'make install' (probably with sudo prefix). This will install Hermes to the target directory of your previous choice (or the default ones).

You can see what is being installed right on the screen (either in bash in Linux, or in Visual Studio output window on Windows).
It is::
    
    the two dynamically-linked libraries ('hermes_common' and 'hermes2d'), these will go to the 'lib' subdirectory of your installation target (and on Windows, the .dll parts to 'bin' subdirectory).
    two include directories ('hermes_common' and 'hermes2d'), these will go to the 'include' subdirectory of your installation target.

The fact that they are linked dynamically is more important on Windows, where you have to **set your PATH environment variable** to point to the directory where you installed it.
If you do not know how to do that, google it.

The recommended way of learning how to use Hermes for your purposes is to take a look at the small collection of what we call 'test examples' right there in the directory.

If you chose not to build them in cmake (using **set(H2D_WITH_TEST_EXAMPLES  NO)**),
change this, and rebuild the library to see them.

If something is not clear enough from the comments in the code, please see the section **Hermes typical example structure**.

These examples (located in hermes2d/test_examples) are::

    00-quickShow
    01-poisson
    02-poisson-newton
    03-navier-stokes
    04-complex-adapt
    05-hcurl-adapt
    06-system-adapt
    07-newton-heat-rk
    08-nonlinearity
    09-trilinos-nonlinear
    10-linear-advection-dg-adapt
    11-transient-adapt
    12-picard
    13-FCT
    14-error-calculation
    15-adaptivity-matrix-reuse-simple
    16-adaptivity-matrix-reuse-layer-interior
    
And these examples are well-documented showcase examples of how to use Hermes.

Take a look at the next section documenting a typical example structure **Hermes typical example structure**.