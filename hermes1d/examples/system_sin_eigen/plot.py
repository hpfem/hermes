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
    for E, eig in eigs:
        if E < -6:
            continue
        if E >= 6:
            break
        x, y = l.get_xy(eig, 0, 100)
        print "plotting E=%f" % E
        plot(x, y)
    show()
