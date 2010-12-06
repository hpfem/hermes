from hermes1d import Linearizer

def plot_file():
    from pylab import plot, show
    import numpy
    data = numpy.loadtxt("solution.gp")
    x = data[:, 0]
    y = data[:, 1]
    plot(x, y)
    show()

def plot_eigs(mesh, eigs, n=10):
    try:
        from jsplot import plot, show
    except ImportError:
        from pylab import plot, show
    import numpy
    l = Linearizer(mesh)
    count = 0
    c = 137.036
    n = 4
    for E, eig in eigs:
        if E < -2*c**2:
            # skip the positron states
            continue
        if E > 0:
            break
        count += 1
        if count > n:
            break
        x, y = l.get_xy(eig, 0, 100)
        print "plotting E=%f" % E
        plot(x, y)
    show()
