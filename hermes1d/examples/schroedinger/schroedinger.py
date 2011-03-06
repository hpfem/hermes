#! /usr/bin/env python

import sys
sys.path.insert(0, "../../../hermes_common")
sys.path.insert(0, "../..")

from numpy import arange, empty, zeros, array
from pylab import plot, show, savefig, grid, gca, legend, figure, title, \
        xlabel, ylabel

from hermes1d import Mesh
from hermes1d.solvers.eigen import solve_eig_numpy, solve_eig_pysparse, \
        solve_eig_scipy
from hermes1d.h1d_wrapper.h1d_wrapper import FESolution, calc_err_est, \
        calc_solution_norm, adapt
from hermes1d.fekete.fekete import Function, Mesh1D
from hermes_common.matrix import CSCMatrix

from _forms import assemble_schroedinger
from plot import plot_eigs, plot_file


def find_element_romanowski(coeffs):
    """
    Finds the smallest coefficient at each element (error) and return the
    element with the largest error.

    Effectivelly it is just using the coefficient at the highest bubble
    function and it needs at least quadratic elements (or higher) to work.
    """
    els = []
    for n, e in enumerate(coeffs):
        error = min(abs(e[2:]))
        #print n, e, error
        els.append((n, error))
    els.sort(key=lambda x: x[1])
    els.reverse()
    n, error = els[0]
    return n, error

def refine_mesh_romanowski(mesh, solutions):
    """
    Uses Romanowski refinement for all solutions in 'solutions'.

    Solutions are given as vectors coming from the matrix solver.
    """
    els2refine = []
    errors = []
    for sol in solutions:
        s = FESolution(mesh, sol)
        id, error = find_element_romanowski(s.get_element_coeffs())
        els2refine.append(id)
        errors.append(error)
    els2refine = list(set(els2refine))
    print "Will refine the elements:", els2refine
    mesh = refine_mesh(mesh, els2refine)
    return mesh

def refine_mesh_h1_adapt(mesh, solutions):
    """
    Uses H1 adaptivity refinement for all solutions in 'solutions'.

    Solutions are given as vectors coming from the matrix solver.
    """
    # so far only for one solution:
    assert len(solutions) == 1
    sol = solutions[0]
    mesh_ref = mesh.reference_refinement()
    print "Fine mesh created (%d DOF)." % mesh_ref.get_n_dof()
    return mesh_ref, [1.0]

def refine_mesh(mesh, els2refine):
    new_pts = []
    pts, orders = mesh.get_mesh_data()
    new_pts.append(pts[0])
    for id in range(len(orders)):
        if id in els2refine:
            new_pts.append((pts[id]+pts[id+1])/2.)
        new_pts.append(pts[id+1])
    # assumes uniform order:
    orders = [orders[0]] * (len(new_pts)-1)
    return Mesh(new_pts, orders)

def do_plot(x, y, n, l):
    n_r = n - l - 1
    styles = {0: "-s", 1: "--o", 2: ":^", 3: "-.v", 4: "-.v", 5: "-.v"}
    plot(x, y, "k" + styles[n_r], label="$R_{%d%d}$" % (n, l))

    grid(True)
    ax = gca()
    xlabel("DOFs")
    ylabel("$E_{num}-E$")
    ax.set_yscale("log")
    title("l=%d" % l)
    legend()

def plot_conv(conv_graph, exact=None, l=None):
    assert exact is not None
    assert l is not None
    n_eig = len(conv_graph[0][1])
    x = []
    y = [[] for n in range(n_eig)]
    for dofs, energies in conv_graph:
        x.append(dofs)
        for i in range(n_eig):
            y[i].append(energies[i]-exact[i])
    f = open("data.py", "w")
    f.write("R_x = {\n")
    f.write("        %d: %s,\n" % (l, x))
    f.write("    }\n")
    f.write("R_y = {\n")
    for i in range(n_eig):
        n = l+1+i
        f.write("        (%d, %d): %s,\n" % (n, l, y[i]))
        #do_plot(x, y[i], n, l)
    f.write("    }\n")
    #savefig("conv_l_0.png")

def flip_vectors(mesh, eigs, mesh_ref, eigs_ref, test_it=False):
    x_c = 1e-3
    for i in range(len(eigs)):
        s = FESolution(mesh, eigs[i])
        s_ref = FESolution(mesh_ref, eigs_ref[i])
        if s.value(x_c) < 0:
            #print "  Multiplying %d-th coarse eigenvector by (-1)" % i
            eigs[i] = -eigs[i]
        if s_ref.value(x_c) < 0:
            #print "  Multiplying %d-th ref. eigenvector by (-1)" % i
            eigs_ref[i] = -eigs_ref[i]

        if test_it:
            # Test it:
            s = FESolution(mesh, eigs[i]).to_discrete_function()
            s_ref = FESolution(mesh_ref, eigs_ref[i]).to_discrete_function()
            same_norm = (s-s_ref).l2_norm()
            flipped_norm = (s+s_ref).l2_norm()
            print same_norm, flipped_norm
            if same_norm > flipped_norm:
                c = min(same_norm, flipped_norm) / max(same_norm, flipped_norm)
                print "Warning: the flip is wrong, c=", c
                # If "c" is almost one, then the vectors can't really be
                # aligned anyway:
                assert c > 0.9
                #s.plot(False)
                #s_ref.plot()

def solve_schroedinger(mesh, l=0, Z=1, eqn_type="R", eig_num=4):
    """
    Solves the Schroedinger equation on the given mesh.

    Returns the energies and eigenfunctions.
    """
    # TODO: return the eigenfunctions as FESolutions
    N_dof = mesh.assign_dofs()
    A = CSCMatrix(N_dof)
    B = CSCMatrix(N_dof)
    assemble_schroedinger(mesh, A, B, l=l, Z=Z, eqn_type=eqn_type)
    eigs = solve_eig_scipy(A.to_scipy_csc(), B.to_scipy_csc())
    eigs = eigs[:eig_num]
    energies = [E for E, eig in eigs]
    #print energies
    assert len(eigs) == eig_num
    eigs = [eig for E, eig in eigs]
    return N_dof, array(energies), eigs

def adapt_mesh(mesh, eigs, l=0, Z=1, adapt_type="hp", eqn_type="R"):
    """
    Adapts the mesh using the adaptivity type 'adapt_type'.

    Returns a new instance of the H1D mesh.

    adapt_type .... one of: h, hp, p, uniform-p, romanowski
    """
    if adapt_type == "romanowski":
        m = refine_mesh_romanowski(mesh, eigs)
        pts, orders = m.get_mesh_data()
        return Mesh(pts, orders)
    elif adapt_type == "uniform-p":
        pts, orders = mesh.get_mesh_data()
        orders = array(orders) + 1
        return Mesh(pts, orders)
    elif adapt_type in ["h", "p", "hp"]:
        NORM = 1 # 1 ... H1; 0 ... L2;
        THRESHOLD = 0.7
        mesh_ref = mesh.reference_refinement()
        print "Fine mesh created (%d DOF)." % mesh_ref.get_n_dof()
        N_dof, energies, eigs_ref = solve_schroedinger(mesh_ref, l=l, Z=Z,
                eqn_type=eqn_type, eig_num=len(eigs))
        flip_vectors(mesh, eigs, mesh_ref, eigs_ref)
        print "    Done."
        sols = []
        sols_ref = []
        print "Normalizing solutions..."
        for i in range(len(eigs)):
            e = (eigs[i]).copy()
            coarse_h1_norm = FESolution(mesh, e).h1_norm()
            e /= coarse_h1_norm
            sols.append(e)
            e = (eigs_ref[i]).copy()
            reference_h1_norm = FESolution(mesh_ref, e).h1_norm()
            e /= reference_h1_norm
            sols_ref.append(e)
            #print "H1 norms:"
            #print "coarse    (%d):" % i, coarse_h1_norm
            #print "reference (%d):" % i, reference_h1_norm
        print "    Done."
        meshes = []
        mesh_orig = mesh.copy()
        mesh_orig.assign_dofs()
        errors = []
        for sol, sol_ref in zip(sols, sols_ref):
            mesh = mesh_orig.copy()
            mesh.assign_dofs()
            mesh_ref = mesh.reference_refinement()
            mesh_ref.assign_dofs()
            mesh.copy_vector_to_mesh(sol, 0)
            mesh_ref.copy_vector_to_mesh(sol_ref, 0)
            err_est_total, err_est_array = calc_error_estimate(NORM, mesh, mesh_ref)
            ref_sol_norm = calc_solution_norm(NORM, mesh_ref)
            err_est_rel = err_est_total/ref_sol_norm
            print "Relative error (est) = %g %%\n" % (100.*err_est_rel)
            errors.append(err_est_rel)
            # TODO: adapt using all the vectors:
            # 0 ... hp, 1 ... h, 2 ... p
            if adapt_type == "hp":
                ADAPT_TYPE = 0
            elif adapt_type == "h":
                ADAPT_TYPE = 1
            elif adapt_type == "p":
                ADAPT_TYPE = 2
            else:
                raise ValueError("Unkown adapt_type")
            adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array, mesh, mesh_ref)
            meshes.append(mesh)
        pts, orders = mesh_orig.get_mesh_data()
        mesh = Mesh1D(pts, orders)
        for m in meshes:
            pts, orders = m.get_mesh_data()
            m = Mesh1D(pts, orders)
            mesh = mesh.union(m)
        pts, orders = mesh.get_mesh_data()
        mesh = Mesh(pts, orders)
        return mesh
    else:
        raise ValueError("Unknown adapt_type")

def create_uniform_mesh(a=0, b=100, n_elem=4):
    """
    Creates a uniform mesh.

    Example::

    >>> create_uniform_mesh(0, 100, 4)
    array([   0.,   25.,   50.,   75.,  100.])
    >>> create_uniform_mesh(0, 100, 5)
    array([   0.,   20.,   40.,   60.,   80.,  100.])
    >>> create_uniform_mesh(100, 200, 4)
    array([ 100.,  125.,  150.,  175.,  200.])

    """
    pts = arange(a, b, float(b-a)/(n_elem))
    pts = list(pts) + [b]
    assert len(pts) == n_elem + 1
    return array(pts)

def create_log_mesh(a=0, b=100, par=20, n_elem=4):
    """
    Creates a logarithmic mesh.

    Example::

    >>> create_log_mesh(0, 100, par=20, n_elem=4)
    array([   0.        ,    3.21724644,   11.95019684,   35.65507127,  100.        ])
    >>> create_log_mesh(0, 100, par=40, n_elem=4)
    array([   0.        ,    1.78202223,    7.87645252,   28.71911092,  100.        ])
    >>> create_log_mesh(0, 100, par=100, n_elem=4)
    array([   0.        ,    0.78625046,    4.43570179,   21.37495437,  100.        ])

    """
    r = par**(1./(n_elem-1))
    pts = [(r**i-1)/(r**n_elem-1)*(b-a)+a for i in range(n_elem+1)]
    assert len(pts) == n_elem + 1
    return array(pts)

def radial_schroedinger_equation_adapt(params, error_tol=1e-8):
    if params["mesh_uniform"]:
        pts = create_uniform_mesh(params["a"], params["b"],
                n_elem=params["el_num"])
    else:
        pts = create_log_mesh(params["a"], params["b"],
                par=params["mesh_par1"],
                n_elem=params["el_num"])
    orders = [params["el_order"]]*(len(pts)-1)
    mesh = Mesh(pts, orders)
    conv_graph = []
    l = params["l"]
    Z = params["Z"]
    N_eig = params["eig_num"]
    exact_energies=[-1.*Z**2/(2*n**2) for n in range(1+l,N_eig+1+l)]
    old_energies = None
    eqn_type = params["eqn_type"]
    try:
        for i in range(10000):
            print "-"*80
            print "adaptivity iteration:", i
            if eqn_type == "rR":
                mesh.set_bc_left_dirichlet(0, 0)
                mesh.set_bc_right_dirichlet(0, 0)
            # Use zero dirichlet for eqn_type="R" as well, just to make sure
            # that we agree with sle1d
            mesh.set_bc_right_dirichlet(0, 0)
            pts, orders = mesh.get_mesh_data()
            print "Current mesh:"
            print pts
            print orders
            #stop
            N_dof, energies, eigs = solve_schroedinger(mesh, l=l, Z=Z,
                    eqn_type=eqn_type, eig_num=N_eig)
            for n in range(1, N_eig+1):
                print "%d  %10.5f" % (n, energies[n-1])
            conv_graph.append((N_dof, energies))
            # This doesn't work well:
            if old_energies is not None:
                err = max(abs(old_energies - energies))
                print "Maximum error in energies:", err
                if err < error_tol:
                    break
            #err = max(abs(energies - exact_energies))
            #print "Maximum error in energies:", err
            #if err < error_tol:
            #    break
            old_energies = energies
        #    exact_energies = array(exact_energies)
        #    print energies - exact_energies
            mesh = adapt_mesh(mesh, eigs, l=l, Z=Z,
                    adapt_type=params["adapt_type"])
    finally:
        pass
        #plot_conv(conv_graph, exact=exact_energies, l=l)
    return [e for e in energies if e < 0]
    #plot_eigs(mesh, zip(energies, eigs))

def calculate_states():
    states = {}
    max_l = 0
    for l in range(100):
        p2 = dict(l=l, Z=47, a=0, b=104.315255921, el_num=6,
                el_order=10, eig_num=6, mesh_uniform=True,
                adapt_type="romanowski", eqn_type="R")
        e = radial_schroedinger_equation_adapt(p2, error_tol=1e-6)
        if e == []:
            break
        states[l] = e
        max_l = l

    print states
    print max_l
    print "Saving to data.py"
    open("data.py", "w").write("""\
max_l = %(max_l)s
states = %(states)s""" % {"max_l": max_l, "states": states})

def main():
    calculate_states()


if __name__ == "__main__":
    main()
