from hermes1d import Linearizer

def plot_file():
    from pylab import plot, show
    import numpy
    data = numpy.loadtxt("solution.gp")
    x = data[:, 0]
    y = data[:, 1]
    plot(x, y)
    show()

def plot_eigs(mesh, eigs):
    #try:
    #    from jsplot import plot, show
    #except ImportError:
    from pylab import plot, show, legend
    import numpy
    from _forms import potential_python
    l = Linearizer(mesh)
    for E, eig in eigs:
        if E >= 0:
            break
        x, y = l.get_xy(eig, 0, 100)
        print "plotting E=%f" % E
        plot(x, y, label="E=%f" % E)
    y = [potential_python(_) for _ in x]
    plot(x, y, "k-", label="V(r)")
    legend()
    show()
