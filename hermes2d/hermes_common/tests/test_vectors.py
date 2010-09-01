from numpy import array, abs

try:
    import hermes_common._hermes_common
    normal_import = True
except ImportError:
    normal_import = False

if normal_import:
    # Running in h1d, h2d or h3d:
    from hermes_common._hermes_common import AVector
else:
    # Running from inside hermes_common
    from _hermes_common import AVector

eps = 1e-10

def _eq(a, b):
    return (abs(a-b) < eps).all()

def test_matrix1():
    m = AVector(5)
    m.add(1, 3.5)
    m.add(2, 4.5)
    m.add(3, 1.5)
    m.add(4, 1.5)
    m.add(2, 1)
    d2 = array([0, 3.5, 5.5, 1.5, 1.5])
    assert (abs(m.to_numpy() - d2) < eps).all()
