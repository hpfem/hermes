#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel, ylim, arange, sin

from sympy.physics.hydrogen import R_nl
from sympy import var


def do_plot(n, l, color="k"):
    x = arange(0, 150, 0.1)
    var("z")
    Z = 47
    f = R_nl(n, l, 1./Z, z)
    y = [f.subs(z, _x) for _x in x]
    plot(x, y, color + "-", lw=2, label="$R_{%d%d}$" % (n, l))

    xlabel("$\\rho$")
    ylabel("$R_{nl}(\\rho)$")
    title("Eigenvectors (l=%d, Z=1)" % l)
    legend()

do_plot(1, 0, color="r")
do_plot(25, 0, color="g")
do_plot(50, 0, color="b")

print "Saving to silver_l_0_vec.png"
savefig("silver_l_0_vec.png")
savefig("silver_l_0_vec.pdf")
