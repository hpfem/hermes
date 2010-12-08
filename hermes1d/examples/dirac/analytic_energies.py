from sympy import var
from sympy.physics.hydrogen import E_nl, E_nl_dirac

def print_exact_energies(N, l, Z):
    """
    Prints lowest N exact energies for the given "l" and "Z".
    """
    print "Z=%d" % Z
    print " n l      Schroed        Dirac (up)     Dirac (down)"
    print "-"*55
    for n in range(l+1, l+N+1):
        print "%2d %d %15.8f %15.8f" % (n, l, E_nl(n, Z=Z).n(),
                E_nl_dirac(n, l, Z=Z)),
        if l > 0:
            print "%15.8f" % E_nl_dirac(n, l, False, Z=Z)
        else:
            print

print_exact_energies(5, 0, 1)
print
print_exact_energies(50, 0, 47)
