# Just import everyting from our cython module:

import os

from _hermes2d import *
from plot import (sln2png, plot_sln_mpl, plot_sln_mayavi, ScalarView, MeshView,
        plot_mesh_mpl)
from runtests import test, doctest
from demos import demo_layer

def get_pxd_include():
    """
    Returns an absolute path to *.pxd files that are needed in order to build
    something against hermes2d.
    """
    this_dir = os.path.abspath(os.path.dirname(__file__))
    include_dir = os.path.join(this_dir, "include")
    return os.path.normpath(include_dir)

def get_include():
    """
    Return the directory in the package that contains the hermes2d/*.h header
    files.

    Extension modules that need to compile against hermes2d should use this
    function to locate the appropriate include directory. Using distutils:

      import hermes2d
      Extension('extension_name', ...
                include_dirs=[hermes2d.get_include()])
    """
    this_dir = os.path.abspath(os.path.dirname(__file__))
    include_dir = os.path.join(this_dir, "..", "..", "..",
            "include", "hermes2d")
    return os.path.normpath(include_dir)

def get_lib():
    """
    Returns the path to the *.so libraries that one needs to link
    against.
    """
    this_dir = os.path.abspath(os.path.dirname(__file__))
    lib_dir = os.path.join(this_dir, "..", "..")
    return os.path.normpath(lib_dir)

def raises(ExpectedException, code):
    """
    Tests, that the "code" raises the ExpectedException exception.

    If so, returns True, otherwise False.
    """
    assert isinstance(code, str)
    import sys
    frame = sys._getframe(1)
    loc = frame.f_locals.copy()
    try:
        exec code in frame.f_globals, loc
    except ExpectedException:
        return True
    return False
