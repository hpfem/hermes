from numpy import array, abs

try:
    import hermes_common._hermes_common
    normal_import = True
except ImportError:
    normal_import = False

if normal_import:
    # Running in h1d, h2d or h3d:
    from hermes_common._hermes_common import CooMatrix, CSRMatrix, CSCMatrix
else:
    # Running from inside hermes_common
    from _hermes_common import CooMatrix, CSRMatrix, CSCMatrix

eps = 1e-10

def _eq(a, b):
    return (abs(a-b) < eps).all()

def _coo_conversions_test(m, d2):
    # test the COO matrix
    d1 = m.to_scipy_coo().todense()
    assert _eq(d1, d2)
    # test conversion from COO
    n = CSRMatrix(m)
    d1 = n.to_scipy_csr().todense()
    assert _eq(d1, d2)
    n = CSCMatrix(m)
    d1 = n.to_scipy_csc().todense()
    assert _eq(d1, d2)
    # test conversion CSC <-> CSR
    n = CSRMatrix(n)
    d1 = n.to_scipy_csr().todense()
    assert _eq(d1, d2)
    n = CSCMatrix(n)
    d1 = n.to_scipy_csc().todense()
    assert _eq(d1, d2)


def test_matrix1():
    m = CooMatrix(5)
    m.add(1, 3, 3.5)
    m.add(2, 3, 4.5)
    m.add(3, 4, 1.5)
    m.add(4, 2, 1.5)
    m.add(2, 3, 1)
    d2 = array([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 3.5, 0],
        [0, 0, 0, 5.5, 0],
        [0, 0, 0, 0, 1.5],
        [0, 0, 1.5, 0, 0],
        ])
    _coo_conversions_test(m, d2)

def test_matrix2():
    m = CooMatrix(5)
    m.add(1, 3, 3.5)
    m.add(2, 3, 4.5)
    m.add(3, 4, 1.5)
    m.add(0, 2, 1.5)
    m.add(2, 3, 1)
    d2 = array([
        [0, 0, 1.5, 0, 0],
        [0, 0, 0, 3.5, 0],
        [0, 0, 0, 5.5, 0],
        [0, 0, 0, 0, 1.5],
        [0, 0, 0, 0, 0],
        ])
    _coo_conversions_test(m, d2)

def test_matrix3():
    m = CooMatrix(5)
    m.add(0, 0, 2)
    m.add(0, 1, 3)
    m.add(1, 0, 3)
    m.add(1, 2, 4)
    m.add(1, 4, 6)
    m.add(2, 1, -1)
    m.add(2, 2, -3)
    m.add(2, 3, 2)
    m.add(3, 2, 1)
    m.add(4, 1, 4)
    m.add(4, 2, 2)
    m.add(4, 4, 1)
    d2 = array([
        [2, 3, 0, 0, 0],
        [3, 0, 4, 0, 6],
        [0,-1,-3, 2, 0],
        [0, 0, 1, 0, 0],
        [0, 4, 2, 0, 1],
        ])
    _coo_conversions_test(m, d2)
    # CSR test:
    m = CSRMatrix(m)
    assert _eq(m.IA, [0, 2, 5, 8, 9, 12])
    assert _eq(m.JA, [0, 1, 0, 2, 4, 1, 2, 3, 2, 1, 2, 4])
    assert _eq(m.A, [2, 3, 3, 4, 6, -1, -3, 2, 1, 4, 2, 1])

    # CSC test:
    m = CSCMatrix(m)
    assert _eq(m.JA, [0, 2, 5, 9, 10, 12])
    assert _eq(m.IA, [0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4])
    assert _eq(m.A, [2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1])

def test_matrix4():
    m = CooMatrix(5, is_complex=True)
    m.add(1, 3, 3.5)
    m.add(2, 3, 4.5)
    m.add(3, 4, 1.5)
    m.add(4, 2, 1.5)
    m.add(2, 3, 1)
    d2 = array([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 3.5, 0],
        [0, 0, 0, 5.5, 0],
        [0, 0, 0, 0, 1.5],
        [0, 0, 1.5, 0, 0],
        ])
    _coo_conversions_test(m, d2)

def test_matrix5():
    m = CooMatrix(5, is_complex=True)
    m.add(1, 3, 3.5)
    m.add(2, 3, 4.5)
    m.add(3, 4, 1.5)
    m.add(0, 2, 1.5)
    m.add(2, 3, 1)
    d2 = array([
        [0, 0, 1.5, 0, 0],
        [0, 0, 0, 3.5, 0],
        [0, 0, 0, 5.5, 0],
        [0, 0, 0, 0, 1.5],
        [0, 0, 0, 0, 0],
        ])
    _coo_conversions_test(m, d2)

def test_matrix6():
    m = CooMatrix(5, is_complex=True)
    m.add(1, 3, 3.5+1j)
    m.add(2, 3, 4.5+2j)
    m.add(3, 4, 1.5+3j)
    m.add(0, 2, 1.5+4j)
    m.add(2, 3, 1)
    d2 = array([
        [0, 0, 1.5+4j, 0, 0],
        [0, 0, 0, 3.5+1j, 0],
        [0, 0, 0, 5.5+2j, 0],
        [0, 0, 0, 0, 1.5+3j],
        [0, 0, 0, 0, 0],
        ])
    _coo_conversions_test(m, d2)

def test_matrix7():
    m = CooMatrix(5, is_complex=True)
    m.add(1, 3, 3.5+1j)
    m.add(2, 3, 4.5+2j)
    m.add(3, 4, 1.5+3j)
    m.add(1, 3, 1.5+4j)
    m.add(2, 3, 1+4j)
    d2 = array([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 5+5j, 0],
        [0, 0, 0, 5.5+6j, 0],
        [0, 0, 0, 0, 1.5+3j],
        [0, 0, 0, 0, 0],
        ])
    _coo_conversions_test(m, d2)
