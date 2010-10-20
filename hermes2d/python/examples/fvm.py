from math import atan2, sin, cos

from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix, lil_matrix
from numpy import array, dot, zeros
from numpy.linalg import norm
import pylab
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from hermes2d import Mesh, set_verbose, MeshView
from hermes2d._numerical_flux import numerical_flux, R as R_const, c_v
from hermes2d.examples import get_GAMM_channel_mesh

marker_bottom = 1
marker_right  = 2
marker_top    = 3
marker_left   = 4


class Plot(object):

    def __init__(self):
        self._num_markers = 5
        self._styles = {
                0: {"color": "black", "lw": 1},
                1: {"color": "blue", "lw": 2},
                2: {"color": "green", "lw": 2},
                3: {"color": "red", "lw": 2},
                4: {"color": "orange", "lw": 2},
            }
        self.color_count = [0]*self._num_markers

    def add_edge(self, p0, p1, normal, marker):
        assert marker < self._num_markers
        self.color_count[marker] += 1
        if self.color_count[marker] == 1:
            pylab.plot([p0[0], p1[0]], [p0[1], p1[1]], label=str(marker),
                    **self._styles[marker])
        else:
            pylab.plot([p0[0], p1[0]], [p0[1], p1[1]], **self._styles[marker])
        # normal
        p0 = (p0+p1)/2
        d = normal * 0.1
        pylab.gca().add_patch(pylab.Arrow(p0[0], p0[1], d[0], d[1],
            width = 0.05, color=self._styles[marker]["color"]))

    def show(self):
        pylab.gca().set_aspect("equal")
        pylab.legend()
        pylab.show()

def plot_state(state_on_elements, mesh):
    patches = []
    colors0 = []
    colors1 = []
    colors2 = []
    colors3 = []
    for e in state_on_elements:
        #print e, state_on_elements[e]
        verts = zeros((4, 2))
        for i in range(4):
            coord = mesh.get_element(e).nodes_vertex[i].coord
            verts[i, 0] = coord[0]
            verts[i, 1] = coord[1]
        polygon = Polygon(verts, True)
        patches.append(polygon)
        colors0.append(state_on_elements[e][0])
        colors1.append(state_on_elements[e][1])
        colors2.append(state_on_elements[e][2])
        colors3.append(state_on_elements[e][3])
    p1 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    p1.set_array(pylab.array(colors0))
    p2 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    p2.set_array(pylab.array(colors1))
    p3 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    p3.set_array(pylab.array(colors2))
    p4 = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
    p4.set_array(pylab.array(colors3))
    fig = pylab.figure()
    ax = fig.add_subplot(221)
    ax.add_collection(p1)
    ax.set_aspect("equal")
    ax.autoscale_view()
    pylab.colorbar(p1)
    ax = fig.add_subplot(222)
    ax.add_collection(p2)
    ax.set_aspect("equal")
    ax.autoscale_view()
    pylab.colorbar(p2)
    ax = fig.add_subplot(223)
    ax.add_collection(p3)
    ax.set_aspect("equal")
    ax.autoscale_view()
    pylab.colorbar(p3)
    ax = fig.add_subplot(224)
    ax.add_collection(p4)
    ax.set_aspect("equal")
    ax.autoscale_view()
    pylab.colorbar(p4)


    pylab.show()

class Edge(object):

    def __init__(self, i, marker, pair_nodes):
        self._pair = pair_nodes
        self._normal = self.calculate_normal()
        self._length = self.calculate_length()
        self._elements = [i]
        self._boundary = True
        self._marker = marker

    def set_second_element(self, i):
        self._elements.append(i)
        self._boundary = False

    def calculate_normal(self):
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        t = array([p1[0]-p0[0], p1[1]-p0[1]])
        t = t/norm(t)
        return array([t[1], -t[0]])

    def calculate_length(self):
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        t = array([p1[0]-p0[0], p1[1]-p0[1]])
        return norm(t)

    @property
    def boundary(self):
        return self._boundary

    @property
    def marker(self):
        return self._marker

    @property
    def normal(self):
        return self._normal

    @property
    def length(self):
        return self._length

    def get_point_0(self):
        return array(self._pair[0].coord)

    def get_point_1(self):
        return array(self._pair[1].coord)

    def get_point_middle(self):
        return (self.get_point_0() + self.get_point_1()) / 2

    @property
    def elements(self):
        return self._elements

    def __str__(self):
        s = "elements: %s, boundary: %s, marker: %s, n=%s" % (self.elements,
                self.boundary, self.marker, self.normal)
        return s

    def __plot__(self, plot):
        p0 = self.get_point_0()
        p1 = self.get_point_1()
        plot.add_edge(p0, p1, self.normal, self.marker)

class Edges(object):

    def __init__(self, mesh):
        self._nodes_dict = mesh.nodes_dict
        self._elements = mesh.elements
        self._nodes = self.extract_nodes(self._elements)
        self._edges = self.extract_edges(mesh)

    def extract_nodes(self, elements):
        nodes = list(set(array(elements).flat))
        nodes.sort()

        _nodes = self._nodes_dict.keys()
        _nodes.sort()
        assert nodes == _nodes
        return nodes

    @property
    def edges(self):
        return self._edges

    def extract_edges(self, mesh):
        edges = {}
        for i in range(mesh.num_elements):
            e = mesh.get_element(i)
            if not e.active:
                continue
            nodes_edge = e.nodes_edge
            nodes_vertex = e.nodes_vertex
            for j, ed in enumerate(nodes_edge):
                jp1 = j + 1
                if jp1 >= len(nodes_edge):
                    jp1 = 0
                pair_nodes = (nodes_vertex[j], nodes_vertex[jp1])
                pair = (pair_nodes[0].id, pair_nodes[1].id)
                if pair in edges:
                    raise Exception("This edge was already processed")
                pair_reversed = (pair[1], pair[0])
                if pair_reversed in edges:
                    edges[pair_reversed].set_second_element(i)
                    continue
                edges[pair] = Edge(i, ed.marker, pair_nodes)
        return edges

    def __str__(self):
        s = "Edges:\n"
        for edge in self._edges:
            s += "%12s:  %s\n" % (edge, self._edges[edge])
        return s

    def __plot__(self, plot):
        for e in self._edges:
            self._edges[e].__plot__(plot)

    def plot(self):
        p = Plot()
        self.__plot__(p)
        p.show()

def T_rot(beta):
    # this is the 3D rotation matrix (2.2.51 in Pavel's master thesis)
    # in 2D without the middle column and row, alpha = 0
    alpha = 0
    return array([
        [1, 0, 0, 0],
        [0, cos(alpha)*cos(beta), sin(beta), 0],
        [0, -cos(alpha)*sin(beta), cos(beta), 0],
        [0, 0, 0, 1]
        ])

def calc_p(w):
    w0, w1, w3, w4 = w
    p = R_const/c_v * (w4 - (w1**2 + w3**2)/(2*w0))
    return p

def calculate_flux(edge, state_on_elements):
    w_l = state_on_elements[edge.elements[0]]
    if edge.boundary:
        if edge.marker in [marker_left, marker_right]:
            w_r = array([1., 50., 0., 1.e5])
        elif edge.marker in [marker_top, marker_bottom]:
            alpha = -atan2(edge.normal[1], edge.normal[0])
            p = calc_p(w_l)
            flux_local = array([0., p, 0., 0.])
            return dot(T_rot(alpha), flux_local)
        else:
            raise Exception("Unhandled boundary")
    else:
        w_r = state_on_elements[edge.elements[1]]
    return numerical_flux(w_l, w_r, edge.normal)

def assembly(edges, state_on_elements, tau):
    elem_contrib = {}
    for e in state_on_elements:
        elem_contrib[e] = zeros((4,))
    for e in edges.edges:
        edge = edges.edges[e]
        flux = -calculate_flux(edge, state_on_elements)
        flux *= edge.length
        #print edge.length, edge, flux, edge.elements
        if edge.boundary:
            elem_contrib[edge.elements[0]] -= flux
        else:
            elem_contrib[edge.elements[0]] -= flux
            elem_contrib[edge.elements[1]] += flux
    dof_map = {}
    for dof, i in enumerate(elem_contrib):
        dof_map[i] = dof
    ndofs = 4*len(dof_map)
    A = lil_matrix((ndofs, ndofs))
    rhs = zeros((ndofs,))
    for e in elem_contrib:
        #print "-"*80
        #print "el_id:", e, "dof:", dof_map[e]
        #print elem_contrib[e]
        for i in range(4):
            C = 1.
            rhs[i*ndofs/4 + dof_map[e]] = C*state_on_elements[e][i]/tau - \
                    elem_contrib[e][i]
            A[i*ndofs/4 + dof_map[e], i*ndofs/4 + dof_map[e]] = C*1./tau
    return A, rhs, dof_map

def set_fvm_solution(x, dof_map):
    ndofs = len(x)
    state_on_elements = {}
    for e in dof_map:
        i = dof_map[e]
        a = zeros((4,))
        for j in range(4):
            a[j] = x[j*ndofs/4 + i]
        state_on_elements[e] = a
    return state_on_elements

def main():
    set_verbose(False)
    mesh = Mesh()
    print "Loading mesh..."
    mesh.load(get_GAMM_channel_mesh())
    #mesh.load("domain-quad.mesh")
    #mesh.refine_element(0, 2)
    mesh.refine_element(1, 2)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()

    print "Constructing edges..."
    nodes = mesh.nodes_dict
    edges = Edges(mesh)
    elements = mesh.elements
    print "Done."

    print "Solving..."
    state_on_elements = {}
    for e in mesh.active_elements:
        state_on_elements[e.id] = array([1., 50., 0., 1.e5])
    #print "initial state"
    #print state_on_elements
    tau = 1e-3
    t = 0.
    for i in range(100):
        A, rhs, dof_map = assembly(edges, state_on_elements, tau)
        #print "A:"
        #print A
        #print "rhs:"
        #print rhs
        #stop
        #print "x:"
        x = spsolve(A, rhs)
        #print x
        #print state_on_elements
        state_on_elements = set_fvm_solution(x, dof_map)
        #print state_on_elements
        t += tau
        print "t = ", t
    plot_state(state_on_elements, mesh)
    #print "state_on_elements:"
    #print state_on_elements
    print "Done."

    #edges.plot()
    #mview = MeshView()
    #mview.show(mesh, lib="mpl", method="orders")

main()
