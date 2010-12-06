from data import max_l, states

orbital_name = ["s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n",
        "o", "p"]

my_states = []
for l in range(max_l):
    n = 1 + l
    for E in states[l]:
        my_states.append((E, n, l))
        n += 1

print "-"*80
total_electrons = 0
for E, n, l in my_states:
    print "(l=%d, n=%d)" % (l, n),
    print "%d%s" % (n, orbital_name[l]),
    print "E=%f" % E,
    deg = 2*(2*l+1)
    print "degeneracy 2*(2*l+1) = %d" % deg,
    total_electrons += deg
    print "total electrons = %d" % total_electrons
    n += 1

print "-"*80
print "Sorted by energy"
my_states.sort(key=lambda state: state[0])
total_electrons = 0
for E, n, l in my_states:
    print "(l=%d, n=%d)" % (l, n),
    print "%d%s" % (n, orbital_name[l]),
    print "E=%f" % E,
    deg = 2*(2*l+1)
    print "degeneracy 2*(2*l+1) = %d" % deg,
    total_electrons += deg
    print "total electrons = %d" % total_electrons
    n += 1
