#! /usr/bin/env python

from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel, ylim, xlim

import hydrogen_uniformpfem
import hydrogen_roman4
import hydrogen_roman6
import hydrogen_roman8
import hydrogen_roman12
import hydrogen_roman17
import hydrogen_hpfem
import hydrogen_pfem
import hydrogen_pfem_uniform_init

def do_plot(x, y, n, l, color="k", label=""):
    n_r = n - l - 1
    if n_r == 0:
        plot(x, y, color + "-o", label=label)
    else:
        plot(x, y, color + "-o")

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$ [Ha]")
    ax.set_yscale("log")
    ylim(ymin=1e-8)
    xlim([0, 90])
    title("Eigenvalues (l=%d, Z=1)" % l)
    legend()

n_eig = 3
l = 0
print "Saving to hydrogen_l_0.png"
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_uniformpfem.R_x[l], hydrogen_uniformpfem.R_y[n, l],
            n, l, "k", "uniform $p$-FEM (L)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_roman4.R_x[l], hydrogen_roman4.R_y[n, l],
            n, l, "m", "$h$-FEM (Romanowski p=4, U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_roman6.R_x[l], hydrogen_roman6.R_y[n, l],
            n, l, "r", "$h$-FEM (Romanowski p=6, U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_roman8.R_x[l], hydrogen_roman8.R_y[n, l],
            n, l, "g", "$h$-FEM (Romanowski p=8, U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_roman12.R_x[l], hydrogen_roman12.R_y[n, l],
            n, l, "b", "$h$-FEM (Romanowski p=12, U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_roman17.R_x[l], hydrogen_roman17.R_y[n, l],
            n, l, "c", "$h$-FEM (Romanowski p=17, U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_pfem.R_x[l], hydrogen_pfem.R_y[n, l],
            n, l, "k", "$p$-FEM (L)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_pfem_uniform_init.R_x[l],
            hydrogen_pfem_uniform_init.R_y[n, l],
            n, l, "k", "$p$-FEM (U)")
for i in range(n_eig):
    n = l+1+i
    do_plot(hydrogen_hpfem.R_x[l], hydrogen_hpfem.R_y[n, l],
            n, l, "k", "$hp$-FEM (U)")
savefig("hydrogen_l_0.png")
savefig("hydrogen_l_0.pdf")
