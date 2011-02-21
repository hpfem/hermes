from _hermes_common cimport c2numpy_double, c2numpy_int

H2D_FN_DX = c_H2D_FN_DX
H2D_FN_DY = c_H2D_FN_DY
H2D_FN_DX_0 = c_H2D_FN_DX_0
H2D_FN_DY_0 = c_H2D_FN_DY_0
H2D_FN_VAL = c_H2D_FN_VAL
H2D_FN_DX = c_H2D_FN_DX
H2D_FN_DY = c_H2D_FN_DY
H2D_FN_DXX = c_H2D_FN_DXX
H2D_FN_DYY = c_H2D_FN_DYY
H2D_FN_DXY = c_H2D_FN_DXY
H2D_FN_DEFAULT = c_H2D_FN_DEFAULT
H2D_FN_ALL = c_H2D_FN_ALL
H2D_EPS_NORMAL = c_H2D_EPS_NORMAL
H2D_EPS_HIGH = c_H2D_EPS_HIGH

H2D_TOTAL_ERROR_REL = c_H2D_TOTAL_ERROR_REL
H2D_TOTAL_ERROR_ABS = c_H2D_TOTAL_ERROR_ABS
H2D_ELEMENT_ERROR_REL = c_H2D_ELEMENT_ERROR_REL
H2D_ELEMENT_ERROR_ABS = c_H2D_ELEMENT_ERROR_ABS

cdef class Nurbs:
    cdef c_Nurbs *thisptr

    @property
    def degree(self):
        return self.thisptr.degree

    @property
    def np(self):
        return self.thisptr.np

    @property
    def pt(self):
        cdef double3 *pt = self.thisptr.pt
        cdef int np = self.thisptr.np
        cdef ndarray vec = c2numpy_double(<double *>pt, 3*np)
        return vec.reshape((np, 3))

    @property
    def nk(self):
        return self.thisptr.nk

    @property
    def kv(self):
        cdef double *kv = self.thisptr.kv
        cdef int nk = self.thisptr.nk
        cdef ndarray vec = c2numpy_double(<double *>kv, nk)
        return vec
        #return vec.reshape((kv,))

    @property
    def ref(self):
        return self.thisptr.nk

cdef class CurvMap:
    cdef c_CurvMap *thisptr

    @property
    def order(self):
        return self.thisptr.order

    @property
    def toplevel(self):
        return self.thisptr.toplevel

    def get_nurbs(self, int k):
        if self.thisptr.toplevel == 0:
            return None
        if self.thisptr.nurbs[k] == NULL:
            return None

        cdef Nurbs n = Nurbs()
        n.thisptr = self.thisptr.nurbs[k]
        return n

cdef class Node:
    cdef c_Node *thisptr

    @property
    def id(self):
        return self.thisptr.id

    @property
    def ref(self):
        return self.thisptr.ref

    @property
    def type(self):
        return self.thisptr.type

    @property
    def bnd(self):
        return self.thisptr.bnd

    @property
    def marker(self):
        return self.thisptr.marker

    @property
    def used(self):
        return self.thisptr.used == 1

    def __str__(self):
        return "Node %d: coord=%r, used=%r" % (self.id, self.coord, self.used)

    @property
    def coord(self):
        return self.thisptr.x, self.thisptr.y

cdef class Element:
    cdef c_Element *thisptr

    @property
    def id(self):
        return self.thisptr.id

    @property
    def nvert(self):
        return self.thisptr.nvert

    @property
    def active(self):
        return self.thisptr.active == 1

    @property
    def used(self):
        return self.thisptr.used == 1

    @property
    def marker(self):
        return self.thisptr.marker

    @property
    def nodes_vertex(self):
        return [self.get_vertex_node(i) for i in range(self.nvert)]

    @property
    def nodes_vertex_id(self):
        return [node.id for node in self.nodes_vertex]

    @property
    def nodes_edge(self):
        return [self.get_edge_node(i) for i in range(self.nvert)]

    def get_vertex_node(self, int id):
        cdef Node n = Node()
        n.thisptr = self.thisptr.vn[id]
        return n

    @property
    def curved_map(self):
        if self.thisptr.cm == NULL:
            return None

        cdef CurvMap cm = CurvMap()
        cm.thisptr = self.thisptr.cm
        return cm

    def get_edge_node(self, int id):
        cdef Node n = Node()
        n.thisptr = self.thisptr.en[id]
        return n

    def get_son_element(self, int id):
        cdef Element e = Element()
        e.thisptr = self.thisptr.sons[id]
        return e

    def __str__(self):
        nodes_id = [node.id for node in self.nodes_vertex]
        return "Element %d: nodes_id=%r, active=%r, marker=%d, used=%r" % \
                (self.id, nodes_id, self.active, self.marker, self.used)

    def get_diameter(self):
        return self.thisptr.get_diameter()

    def get_area(self):
        return self.thisptr.get_area()

def get_node_id(node):
    # This function is used in nodes.sort() in the Mesh class.
    # Cython doesn't yet support lambda functions, nor closures.
    return node.id


cdef class Mesh:

    def __cinit__(self):
        self.thisptr = new_Mesh()

    def __dealloc__(self):
        delete(self.thisptr)

    def copy(self, Mesh m):
        self.thisptr.copy(m.thisptr)

    def load(self, char* filename):
        cdef c_H2DReader mloader
        mloader.load(filename, self.thisptr)

    def load_str(self, char* mesh):
        cdef c_H2DReader mloader
        mloader.load_str(mesh, self.thisptr)

    @property
    def elements_markers(self):
        """
        Returns the list of elements (as ids) made up by their corresponding vertices,
        plus an extra number denoting the element (material) marker.

        Element markers allow you to use different material parameters in areas 
        with different material identity.
 
        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh() 
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.elements_markers
        [[0, 1, 4, 3, 0], [3, 4, 7, 0], [3, 7, 6, 0], [2, 3, 6, 5, 0]] 

        If the domain is composed of only one material, as in the example above,
        all elements may be assigned a zero marker.  In the example above
        "[0, 1, 4, 3, 0]" is an element with "0" being the corresponding element
        marker.
    
        """
        element_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            element_list.append(el.nodes_vertex_id + [el.marker])
        return element_list

    @property
    def elements(self):
        """
        Returns the list of elements (as ids) made up by their corresponding vertices.

        Elements in a mesh are made up of zero based indices of their vertices in a 
        counterclockwise order.

        Example:
        >>> import hermes2d
        >>> m = hermes2d.Mesh() 
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.elements
        [[0, 1, 4, 3], [3, 4, 7], [3, 7, 6], [2, 3, 6, 5]] 

        In the example above "[0, 1, 4, 3]" is an element.

        """
        element_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                element_list.append(el.nodes_vertex_id)
        return element_list

    @property
    def active_elements(self):
        """
        Return the list of active elements (as Element instances).
        """
        element_list = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                element_list.append(el)
        return element_list

    @property
    def nodes(self):
        """
        Returns the list of nodes made up by the coordinates of the vertices.

        Note that this function works even after mesh refinement.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh() 
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.nodes
        [(0.0, -1.0), (1.0, -1.0), (-1.0, 0.0), (0.0, 0.0), (1.0, 0.0), (-1.0,
        1.0), (0.0, 1.0), (0.70710700000000004, 0.70710700000000004)]

        >>> m.refine_all_elements(); 
        >>> m.nodes
        [(0.0, -1.0), (1.0, -1.0), (-1.0, 0.0), (0.0, 0.0), (1.0, 0.0), (-1.0,
	1.0), (0.0, 1.0), (0.70710700000000004, 0.70710700000000004), (0.5,
	0.0), (1.0, -0.5), (0.0, 0.5), (0.5, -1.0), (0.38268352000946509,
	0.92387966368036412), (-1.0, 0.5), (-0.5, 0.5), (-0.5, 1.0), (-0.5,
	0.0), (0.0, -0.5), (0.5, -0.5), (0.92387966368036412,
	0.38268352000946509), (0.35355350000000002, 0.35355350000000002)]

        Notice how in the example above we refined our mesh with "m.refine_all_elements();"
        and diplayed our new "list of nodes" after refinement.  In the example above
        "(0.0, -1.0)" represents a node.       

        """
        # This is a really slow implementation, but it gets the job
        # done for now. Later on, we should get the list of nodes from C++
        # directly.
        nodes = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            nodes.extend(el.nodes_vertex)
        node_dict = {}
        for node in nodes:
            node_dict[node.id] = node
        nodes = node_dict.values()
        nodes.sort(key=get_node_id)
        nodes_coord = [node.coord for node in nodes]
        return nodes_coord

    @property
    def nodes_dict(self):
        """
        Returns a Python dictionary of all the active nodes made up by the
        coordinates of their corresponding vertices.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh() 
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.nodes_dict
        {0: (0.0, -1.0), 1: (1.0, -1.0), 2: (-1.0, 0.0), 3: (0.0, 0.0), 4: (1.0,
        0.0), 5: (-1.0, 1.0), 6: (0.0, 1.0), 7: (0.70710700000000004,
        0.70710700000000004)}

        In the example above "0: (0.0, -1.0)", "(0.0, -1.0)" is a node and the lone "0"
        to the left of the colon is the corresponding dictionary item number.

        """
        nodes = []
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                nodes.extend(el.nodes_vertex)
        node_dict = {}
        for node in nodes:
            node_dict[node.id] = node.coord
        return node_dict

    @property
    def curves(self):
        """
        Return curves
        """
        cdef Element e 
        cdef double2 *tp
        cdef ndarray vec

        crv = {}
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                cm = el.curved_map
                if cm != None:
                    if cm.toplevel == True:
                        for j in range(4):
                            n = cm.get_nurbs(j)
                            if n != None:
                                sp = (n.pt[0][0], n.pt[0][1])
                                cp = (n.pt[1][0], n.pt[1][1])
                                ep = (n.pt[2][0], n.pt[2][1])
                                crv[i] = []
                                crv[i].append([sp, cp, ep])
                    else:
                        e = el
                        tp = transform(e.thisptr)
                        vec = c2numpy_double(<double *>tp, 2*5)
                        crv[i] = []
                        x = []
                        for point in vec.reshape((5, 2)):
                            x.append((point[0], point[1]))
                        crv[i].append(x)
        print crv
        return crv

    def get_polygonal_boundary(self):
        """
        Return the polygonal approximation of boundaries of all elements.
        """
        cdef Element e
        cdef double2 *tp
        cdef int npoints
        cdef ndarray vec

        crv = {}
        for i in range(self.num_elements):
            e = self.get_element(i)
            if e.active:
                element_polygonal_boundary(e.thisptr, &tp, &npoints)
                vec = c2numpy_double(<double *>tp, 2*npoints)
                crv[e.id] = vec.reshape((npoints, 2))
        return crv

    @property
    def num_elements(self):
        return self.thisptr.get_num_elements()

    @property
    def num_base_elements(self):
        return self.thisptr.get_num_base_elements()

    @property
    def num_active_elements(self):
        return self.thisptr.get_num_active_elements()

    @property
    def max_element_id(self):
        return self.thisptr.get_max_element_id()

    def __str__(self):
        return "Mesh with %d elements: active=%d, base=%d, max_id=%d" % \
                (self.num_elements, self.num_active_elements,
                        self.num_base_elements, self.max_element_id)

    def load_hermes(self, char* filename):
        """
        Loads a mesh (in the hermes format) from a file.

        This uses a pure Python reader. If you want to use the C++ flex reader,
        use load().
        """
        from mesh import read_hermes_format
        nodes, elements, boundary, nurbs = \
                read_hermes_format(filename)
        self.create(nodes, elements, boundary, nurbs)

    def get_elements_order(self, space):
        """
        Returns list of orders
        """
        orders_list = {}
        for i in range(self.num_elements):
            el = self.get_element(i)
            if el.active:
                order = space.get_element_order(i)
                h = order & ((1 << 5) - 1)
                v = order >> 5

                import math
                # This uses the maximum of "h", "v", as we can't yet plot
                # anisotropic polynomial degrees
                ord = max(h, v)
                #ord = int(((h+v)/2.0))
                if ord == 0:
                    ord = 1

                #orders_list.append(int(((h+v)/2.0)))
                orders_list[el.id] = ord

        return orders_list

    def create(self, nodes, elements, boundary, nurbs):
        """
        Creates a mesh from a list of nodes, elements, boundaries and nurbs.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh()
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])

        """
        if boundary is None:
            boundary = []
        if nurbs is None:
            nurbs = []
        m = "1 0\n"
        m += "%d\n" % len(nodes)
        for node in nodes:
            m += "%f %f\n" % tuple(node)
        m += "%d\n" % len(elements)
        for el in elements:
            m += ("%d "*len(el)) % tuple(el)
            m += "\n"
        m += "%d\n" % len(boundary)
        for b in boundary:
            m += "%d %d %d\n" % tuple(b)
        m += "%d\n" % len(nurbs)
        for n in nurbs:
            m += "%d %d\n0\n%f\n" % tuple(n)
        m += "0\n"
        self.load_str(m)

    def save(self, char* filename):
        self.thisptr.save(filename)

    def refine_element_id(self, int id, int refinement=0):
        """
        Refines the element of interest.

        The first parameter takes your element of interest, and the second 
        parameter takes your refinement option.  0:isotropically  1:anisotropically

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh()
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.refine_element_id(0,0);
        >>> m.elements 
        [[3, 4, 7], [3, 7, 6], [2, 3, 6, 5], [0, 11, 20, 19], [11, 1, 9, 20],
        [20, 9, 4, 8], [19, 20, 8, 3]]

        As can be seen from the example above, more elements were created after the
        refinement of the element of interest; in this case element "0" with refinement
        option "0".  Originally we had four elements, but after the refinement of the
        element of interest we now have seven.               

        """
        self.thisptr.refine_element_id(id, refinement)

    def refine_all_elements(self):
        """
        Refines all initial elements in the mesh.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh()
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.elements
        [[0, 1, 4, 3], [3, 4, 7], [3, 7, 6], [2, 3, 6, 5]]
        >>> m.refine_all_elements();
        >>> m.elements
        [[0, 11, 20, 19], [11, 1, 9, 20], [20, 9, 4, 8], [19, 20, 8, 3], [3, 8,
        34], [8, 4, 33], [34, 33, 7], [33, 34, 8], [3, 34, 10], [34, 7, 12],
        [10, 12, 6], [12, 10, 34], [2, 18, 16, 15], [18, 3, 10, 16], [16, 10, 6,
        17], [15, 16, 17, 5]] 

        Notice in the example above how we showed the initial elements that were in the mesh
        with the command "m.elements".  We then went ahead and refined these initial elements
        with "m.refine_all_elements();" and displayed the newly refined set of elements with 
        "m.elements".  In the example above, as an example, "[0, 11, 20, 19]" is an element
        made up by its corresponding vertices.

        """
        self.thisptr.refine_all_elements()

    def refine_towards_boundary(self, int marker, int depth):
        self.thisptr.refine_towards_boundary(marker, depth)

    def refine_towards_vertex(self, int marker, int depth):
        """
        Refines a mesh towards a given vertex by a certain refinement degree.

        The first parameter inserted is the vertex of interest, followed by a number
        representing the wanted number of mesh refinements.

        Example:

        >>> import hermes2d
        >>> m = hermes2d.Mesh()
        >>> m.create([
        ...         [0, -1],
        ...         [1, -1],
        ...         [-1, 0],
        ...         [0, 0],
        ...         [1, 0],
        ...         [-1, 1],
        ...         [0, 1],
        ...         [0.707106781, 0.707106781],
        ...     ], [
        ...         [0, 1, 4, 3, 0],
        ...         [3, 4, 7, 0],
        ...         [3, 7, 6, 0],
        ...         [2, 3, 6, 5, 0],
        ...     ], [
        ...         [0, 1, 1],
        ...         [1, 4, 2],
        ...         [3, 0, 4],
        ...         [4, 7, 2],
        ...         [7, 6, 2],
        ...         [2, 3, 4],
        ...         [6, 5, 2],
        ...         [5, 2, 3],
        ...     ], [
        ...         [4, 7, 45],
        ...         [7, 6, 45],
        ...     ])
        >>> m.elements
        [[0, 1, 4, 3], [3, 4, 7], [3, 7, 6], [2, 3, 6, 5]]
        >>> m.refine_towards_vertex(0,2);
        >>> m.elements
        [[3, 4, 7], [3, 7, 6], [2, 3, 6, 5], [11, 1, 9, 20], [20, 9, 4, 8], [19,
        20, 8, 3], [0, 24, 35, 34], [24, 11, 21, 35], [35, 21, 20, 33], [34, 35,
        33, 19]]
        
        In the example above, to illustrate the post mesh refinement towards the given
        vertex of interest,  we first show the  initial mesh elements with the command
        "m.elements".  We then refine our mesh towards vertex "0" (coordinates [0,-1])
        with a refinement multiplication of "2".  This is shown in the command
        "m.refine_towards_vertex(0,2);".  And finally the new list of  elements after
        refinement is shown with "m.elements".  In the example above, as an example,
        "[3, 4, 7]" is an  element made up of its corresponding vertices.  

        """
        self.thisptr.refine_towards_vertex(marker, depth)

    def get_element(self, int id):
        cdef Element e = Element()
        e.thisptr = self.thisptr.get_element(id)
        return e

    def plot(self, *args, **kwargs):
        """
        Plots the mesh and shows it to the user.

        It passes all arguments to the MeshView.show() function, so read its
        documentation for the meaning.
        """
        from hermes2d import MeshView
        m = MeshView()
        m.show(self, *args, **kwargs)

    def show(self, *args, **kwargs):
        """
        Plots the mesh and shows it to the user.

        It passes all arguments to the MeshView.show() function, so read its
        documentation for the meaning.
        """
        self.plot(*args, **kwargs)

    def convert_triangles_to_quads(self):
        self.thisptr.convert_triangles_to_quads()

    def convert_quads_to_triangles(self):
        self.thisptr.convert_quads_to_triangles()

cdef class H1Shapeset:
    cdef c_H1Shapeset *thisptr

    def __cinit__(self):
        self.thisptr = new_H1Shapeset()

    def __dealloc__(self):
        delete(self.thisptr)

cdef class L2Shapeset:
    """ L2 Shapeset.

        Suggested Use: ```l2_shapeset = L2Shapeset()```
    """
    cdef c_L2Shapeset *thisptr

    def __cinit__(self):
        self.thisptr = new_L2Shapeset()

    def __dealloc__(self):
        delete(self.thisptr)

cdef class PrecalcShapeset:
    cdef c_PrecalcShapeset *thisptr

    def __cinit__(self, H1Shapeset h):
        self.thisptr = new_PrecalcShapeset(h.thisptr)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class H1Space:

    def __init__(self, Mesh m = None, p_init = 1, H1Shapeset h = None):
        self.thisptr = new_H1Space(m.thisptr, NULL, NULL, p_init, NULL)

    def __dealloc__(self):
        delete(self.thisptr)

    def set_uniform_order(self, int tri_order):
        self.thisptr.set_uniform_order(tri_order)

    def assign_dofs(self, first_dof=0, stride=1):
        return self.thisptr.assign_dofs(first_dof, stride)

    def copy_orders(self, H1Space s, int inc=0):
        self.thisptr.copy_orders(s.thisptr, inc)

    def get_element_order(self, int el_id):
        return self.thisptr.get_element_order(el_id)

    def get_num_dofs(self):
        return self.thisptr.get_num_dofs()

cdef class L2Space:
    """ L2 Space.

        Suggested Use: ```l2_space = L2Space(mesh, l2_shapeset)```
    """
    def __init__(self, Mesh m):
        self.thisptr = new_L2Space(m.thisptr)

    def __dealloc__(self):
        delete(self.thisptr)

    def set_uniform_order(self, int tri_order):
        """ Sets an order of all elements. In a case of quads, this supplied
        order is set to both the horizontal and the vertical direction

        Suggested Use: ```l2_space.set_uniform_order(3)```
    """
        self.thisptr.set_uniform_order(tri_order)

    def assign_dofs(self, first_dof=0, stride=1):
        """ Assign DOFs. Updates internal structures after orders has been changed
        Returns a number of DOFs.

        Suggested Use: ```l2_space.assign_dofs()```
    """
        return self.thisptr.assign_dofs(first_dof, stride)

    def copy_orders(self, L2Space s, int inc=0):
        """ Copies order from a supplied L2 space.

        Suggested Use: ```l2_space.copy_orders(another_l2_space)```
    """
        self.thisptr.copy_orders(s.thisptr, inc)

    def get_element_order(self, int el_id):
        """ Returns order of element. In case of quadrilaterals, the order
        is in an encoded form.

        Suggested Use: ```order = l2_space.get_element_order(element_id)```
    """
        return self.thisptr.get_element_order(el_id)

    def get_num_dofs(self):
        """ Returns a number of DOFs.

        Suggested Use: ```num_dofs = l2_space.get_num_dofs()```
    """
        return self.thisptr.get_num_dofs()

cdef api object H1Space_from_C(c_H1Space *h):
    cdef H1Space n
    n = <H1Space>PY_NEW(H1Space)
    n.thisptr = h
    return n

cdef api object Mesh_from_C(c_Mesh *h):
    cdef Mesh n
    n = <Mesh>PY_NEW(Mesh)
    n.thisptr = h
    return n

cdef api object LinSystem_from_C(c_LinSystem *h):
    cdef LinSystem n
    n = <LinSystem>PY_NEW(LinSystem)
    n.thisptr = h
    return n

cdef class Transformable:
    pass

cdef class Function(Transformable):
    pass

cdef class ScalarFunction(Function):
    pass

cdef class MeshFunction(ScalarFunction):

    def __add__(x, y):
        return SumFilter(x, y)

    def __sub__(x, y):
        return DiffFilter(x, y)

    def __pow__(x, y, z):
        if y == 2:
            return SquareFilter(x)
        return NotImplemented

    def __neg__(x):
        return x-x-x

    #def integrate(self):
    #    cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
    #    return integrate(m)

    def l2_norm(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return l2_norm(m)

    def h1_norm(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return h1_norm(m)

    def get_pt_value(self, x, y):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return m.get_pt_value(x, y)

    def get_mesh(self):
        cdef c_MeshFunction *m = <c_MeshFunction *>(self.thisptr)
        return Mesh_from_C(m.get_mesh())

cdef api object Solution_from_C(c_Solution *s):
    cdef Solution n
    n = <Solution>PY_NEW(Solution)
    n.thisptr = <c_Function *>s
    return n

cdef class Solution(MeshFunction):

    def __cinit__(self):
        self.thisptr = <c_Function *>new_Solution()

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def copy(self, Solution s):
        (<c_Solution *>(self.thisptr)).copy((<c_Solution *>(s.thisptr)))

    def copy_by_reference(self, Solution s):
        self.thisptr = s.thisptr

    def set_zero(self, Mesh m):
        (<c_Solution *>(self.thisptr)).set_zero(m.thisptr)

    def set_const(self, Mesh m, scalar):
        (<c_Solution *>(self.thisptr)).set_const(m.thisptr, scalar)


    def set_fe_solution(self, H1Space s, PrecalcShapeset pss, ndarray v):
        """
        Sets the solution using the coefficient vector Y.
        """
        cdef int n = len(v)
        cdef int i
        from numpy import array
        cdef ndarray vec = array(v, dtype="double")
        cdef scalar *pvec = <scalar *>vec.data
        (<c_Solution *>(self.thisptr)).set_fe_solution(s.thisptr, pss.thisptr,
                pvec)

    def plot(self, *args, **kwargs):
        """
        Plots the solution and shows it to the user.

        It passes all arguments to the ScalarView.show() function, so read its
        documentation for the meaning.
        """
        from hermes2d import ScalarView
        sview = ScalarView()
        sview.show(self, *args, **kwargs)

    def show(self, *args, **kwargs):
        """
        Plots the solution and shows it to the user.

        It passes all arguments to the ScalarView.show() function, so read its
        documentation for the meaning.
        """
        self.plot(*args, **kwargs)
 
    # the get_fe_solution() method is is not yet implemented in the C++ hermes:
    #def get_fe_solution(self):
    #    """
    #    Returns the Y coefficients vector as a numpy array.
    #    """
    #    cdef int Ylen
    #    cdef scalar *Y
    #    (<c_Solution *>(self.thisptr)).get_fe_solution(&Ylen, &Y)
    #    if not (Ylen > 0 and Y != NULL):
    #        raise Exception("Ylen (%d) or Y is not valid." % Ylen)
    #    from numpy import empty
    #    cdef ndarray vec = empty([Ylen], dtype="double")
    #    cdef scalar *pvec = <scalar *>vec.data
    #    memcpy(pvec, Y, Ylen*sizeof(scalar))
    #    return vec

cdef class Filter(MeshFunction):
    pass

cdef class SimpleFilter(Filter):
    pass

cdef class VonMisesFilter(Filter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, double l, double m):
        self.thisptr = <c_Function *>new_VonMisesFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, l, m)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class MagFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2, int item1, int item2):
        self.thisptr = <c_Function *>new_MagFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class DiffFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=H2D_FN_VAL, int item2=H2D_FN_VAL):
        self.thisptr = <c_Function *>new_DiffFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class SumFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln1, MeshFunction sln2,
            int item1=H2D_FN_VAL, int item2=H2D_FN_VAL):
        self.thisptr = <c_Function *>new_SumFilter(<c_MeshFunction *>sln1.thisptr, <c_MeshFunction *>sln2.thisptr, item1, item2)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class SquareFilter(SimpleFilter):

    def __cinit__(self, MeshFunction sln, int item=H2D_FN_VAL):
        self.thisptr = <c_Function *>new_SquareFilter(<c_MeshFunction *>sln.thisptr, item)

    #def __dealloc__(self):
    #    delete(self.thisptr)

cdef class WeakForm:

    def __cinit__(self, int neq=1):
        self.thisptr = new_WeakForm(neq)

    def __dealloc__(self):
        delete(self.thisptr)

cdef class DummySolver(CommonSolver):

    def __cinit__(self):
        self.thisptr = <c_CommonSolver *>(new_DummySolver())

    def __dealloc__(self):
        delete(self.thisptr)

cdef class LinSystem:

    def __init__(self, WeakForm wf):
        self.thisptr = new_LinSystem(wf.thisptr)
   # def __init__(self, WeakForm wf, CommonSolver solver):
   #     self.thisptr = new_LinSystem(wf.thisptr, solver.thisptr)
   # def __init__(self, WeakForm wf, H1Space sp):
   #    self.thisptr = new_LinSystem(wf.thisptr, sp.thisptr)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def set_spaces(self, *args):
        self._spaces = args
        cdef int n = len(args)
        cdef H1Space a, b, c, d
        if n == 1:
            a = args[0]
            self.thisptr.set_spaces2(n, a.thisptr)
        elif n == 2:
            a, b = args
            self.thisptr.set_spaces2(n, a.thisptr, b.thisptr)
        elif n == 3:
            a, b, c = args
            self.thisptr.set_spaces2(n, a.thisptr, b.thisptr, c.thisptr)
        elif n == 4:
            a, b, c, d = args
            self.thisptr.set_spaces2(n, a.thisptr, b.thisptr, c.thisptr,
                    d.thisptr)
        else:
            raise NotImplementedError()

    def set_pss(self, *args):
        """
        self._pss = args
        cdef int n = len(args)
        cdef PrecalcShapeset s1, s2, s3, s4
        if n == 1:
            s1, = args
            self.thisptr.set_pss(n, s1.thisptr)
        elif n == 2:
            s1, s2 = args
            self.thisptr.set_pss(n, s1.thisptr, s2.thisptr)
        elif n == 3:
            s1, s2, s3 = args
            self.thisptr.set_pss(n, s1.thisptr, s2.thisptr, s3.thisptr)
        elif n == 4:
            s1, s2, s3, s4 = args
            self.thisptr.set_pss(n, s1.thisptr, s2.thisptr, s3.thisptr,
                    s4.thisptr)
        else:
            raise NotImplementedError()
        """
        pass
        
    def set_pss2(self, *args):
        self._pss = args
        cdef int n = len(args)
        cdef PrecalcShapeset s1, s2, s3, s4
        if n == 1:
            s1, = args
            self.thisptr.set_pss2(n, s1.thisptr)
        elif n == 2:
            s1, s2 = args
            self.thisptr.set_pss2(n, s1.thisptr, s2.thisptr)
        elif n == 3:
            s1, s2, s3 = args
            self.thisptr.set_pss2(n, s1.thisptr, s2.thisptr, s3.thisptr)
        elif n == 4:
            s1, s2, s3, s4 = args
            self.thisptr.set_pss2(n, s1.thisptr, s2.thisptr, s3.thisptr,
                    s4.thisptr)
        else:
            raise NotImplementedError()


    def solve_system(self, *args, lib="scipy"):
        """
        Solves the linear system using scipy.
        """
        cdef int n = len(args)

        cdef Solution s0, s1, s2, s3
        cdef ndarray vec
        cdef scalar *pvec
        
        """
        if lib == "hermes":

            if n == 1:
                s0 = args[0]
                self.thisptr.solve(n, s0.thisptr)
            elif n == 2:
                s0, s1 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr)
            elif n == 3:
                s0, s1, s2 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr, s2.thisptr)
            elif n == 4:
                s0, s1, s2, s3 = args
                self.thisptr.solve(n, s0.thisptr, s1.thisptr, s2.thisptr,
                        s3.thisptr)
            else:
                raise NotImplementedError()
        """
        
        if lib == "hermes":

            if n == 1:
                s0 = args[0]
                self.thisptr.solve2(n, s0.thisptr)
            elif n == 2:
                s0, s1 = args
                self.thisptr.solve2(n, s0.thisptr, s1.thisptr)
            elif n == 3:
                s0, s1, s2 = args
                self.thisptr.solve2(n, s0.thisptr, s1.thisptr, s2.thisptr)
            elif n == 4:
                s0, s1, s2, s3 = args
                self.thisptr.solve2(n, s0.thisptr, s1.thisptr, s2.thisptr,
                        s3.thisptr)
            else:
                raise NotImplementedError()
                          
        elif lib == "scipy":
            import warnings
            from scipy.sparse.linalg import cg
            from scipy.sparse.linalg import spsolve
            from numpy import array
            A = self.get_matrix()
            rhs = self.get_rhs()
            #x, res = cg(A, rhs)
            # Ignore the warning in scipy:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                x = spsolve(A, rhs)
            vec = array(x, dtype="double")
            pvec = <scalar *>vec.data

            for i, sln in enumerate(args):
                (<c_Solution *>((<Solution>sln).thisptr)).set_fe_solution(
                        self.thisptr.get_space(i),
                        self.thisptr.get_pss(i),
                        pvec)
        else:
            raise NotImplementedError("Unknown library")

    def assemble(self):
        self.thisptr.assemble()

    def get_matrix_csc(self):
        """
        Returns the matrix A as a (Ap, Ai, Ax) tuple.

        See also DiscreteProblem.get_matrix() to get a SciPy matrix.
        """
        cdef int *Ap, *Ai, n, nnz
        cdef scalar *Ax
        self.thisptr.get_matrix(Ap, Ai, Ax, n)
        nnz = Ap[n]
        aAp = c2numpy_int(Ap, n+1)
        aAi = c2numpy_int(Ai, nnz)
        aAx = c2numpy_double(Ax, nnz)
        return aAp, aAi, aAx

    def get_matrix(self):
        """
        Returns the global matrix A as a SciPy matrix.
        """
        from scipy.sparse import csc_matrix
        Ap, Ai, Ax = self.get_matrix_csc()
        return csc_matrix((Ax, Ai, Ap))

    def get_rhs(self):
        """
        Return the RHS as a numpy array
        """
        cdef scalar *rhs
        cdef int n
        self.thisptr.get_rhs(rhs, n)
        return c2numpy_double(rhs, n)

    def get_num_dofs(self):
        self.thisptr.get_num_dofs()

cdef class RefSystem(LinSystem):

    def __init__(self, LinSystem ls):
        self.thisptr = <c_LinSystem *>new_RefSystem(ls.thisptr)

    def assemble(self):
        (<c_RefSystem *>(self.thisptr)).assemble()

    # this is commented out, because get_ref_space() is not yet implemented in
    # C++ hermes2d
    #def get_ref_space(self, int eq):
    #    cdef c_H1Space *r = <c_H1Space *>(
    #            (<c_RefSystem *>(self.thisptr)).get_ref_space(eq)
    #        )
    #    return H1Space_from_C(r)


#cdef class DiscreteProblem:
#
#    def __cinit__(self):
#        self.thisptr = new_DiscreteProblem()
#
#    def __dealloc__(self):
#        delete(self.thisptr)
#
#    def copy(self, DiscreteProblem ep):
#        self.thisptr.copy(ep.thisptr)
#
#    def set_num_equations(self, int neq):
#        self.thisptr.set_num_equations(neq)
#
#    def set_spaces(self, *args):
#        cdef int n = len(args)
#        cdef H1Space a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.set_spaces(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.set_spaces(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.set_spaces(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def set_external_fns(self, *args):
#        cdef int n = len(args)
#        cdef Solution a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.set_external_fns(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.set_external_fns(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def set_pss(self, *args):
#        cdef int n = len(args)
#        cdef PrecalcShapeset s
#        if n == 1:
#            s = args[0]
#            self.thisptr.set_pss(n, s.thisptr)
#        else:
#            raise NotImplementedError()
#
#    #def set_bilinear_form(self, int i, int j, BiFormFnVol unsym):
#    #    self.thisptr.set_bilinear_form(i, j, unsym)
#
#    def create_matrix(self):
#        self.thisptr.create_matrix()
#
#    def set_quiet(self, quiet=True):
#        self.thisptr.set_quiet(quiet)
#
#    def assemble_matrix_and_rhs(self):
#        self.thisptr.assemble_matrix_and_rhs()
#
#    def solve_system(self, *args):
#        cdef int n = len(args)
#        cdef Solution a, b, c
#        if n == 1:
#            a = args[0]
#            self.thisptr.solve_system(n, a.thisptr)
#        elif n == 2:
#            a, b = args
#            self.thisptr.solve_system(n, a.thisptr, b.thisptr)
#        elif n == 3:
#            a, b, c = args
#            self.thisptr.solve_system(n, a.thisptr, b.thisptr, c.thisptr)
#        else:
#            raise NotImplementedError()
#
#    def save_matrix_matlab(self, char *filename, char *varname):
#        self.thisptr.save_matrix_matlab(filename, varname)
#
#    def save_matrix_coo(self, char *filename):
#        self.thisptr.save_matrix_coo(filename)
#
#    def get_matrix_csc(self):
#        """
#        Returns the matrix A as a (Ap, Ai, Ax) tuple.
#
#        See also DiscreteProblem.get_matrix() to get a SciPy matrix.
#        """
#        cdef int *Ap, *Ai, n, nnz
#        cdef scalar *Ax
#        self.thisptr.get_matrix(Ap, Ai, Ax, n)
#        nnz = Ap[n]
#        from numpy import zeros
#        cdef ndarray aAp = zeros(n+1, dtype="int32")
#        cdef ndarray aAi = zeros(nnz, dtype="int32")
#        cdef ndarray aAx = zeros(nnz, dtype="double")
#        cdef int *pAp = <int *>aAp.data
#        cdef int *pAi = <int *>aAi.data
#        cdef double *pAx = <double *>aAx.data
#        cdef int i
#        for i in range(n+1):
#            pAp[i] = Ap[i]
#        for i in range(nnz):
#            pAi[i] = Ai[i]
#            pAx[i] = Ax[i]
#        return aAp, aAi, aAx
#
#    def get_matrix(self):
#        """
#        Returns the global matrix A as a SciPy matrix.
#        """
#        from scipy.sparse import csc_matrix
#        Ap, Ai, Ax = self.get_matrix_csc()
#        return csc_matrix((Ax, Ai, Ap))
#
#    def free_matrix(self):
#        self.thisptr.free_matrix()

cdef class CandList:
    """ A prefined list of andidates.

        Possible values are attributes of the class:

        - P_ISO: P-candidates only. Orders are modified uniformly.
    - P_ANISO: P-candidates only. Orders are modified non-uniformly.
    - H_ISO: H-candidates only. Orders are not modified.
    - H_ANISO: H- and ANISO-candidates only. Orders are not modified.
    - HP_ISO: H- and P-candidates only. Orders are modified uniformly.
    - HP_ANISO_H: H-, ANISO- and P-candidates. Orders are modified uniformly.
    - HP_ANISO_P: H- and P-candidates only. Orders are modified non-uniformly.
    - HP_ANISO: H-, ANISO- and P-candidates. Orders are modified non-uniformly.
    """
    H2D_P_ISO = c_H2D_P_ISO
    H2D_P_ANISO = c_H2D_P_ANISO
    H2D_H_ISO = c_H2D_H_ISO
    H2D_H_ANISO = c_H2D_H_ANISO
    H2D_HP_ISO = c_H2D_HP_ISO
    H2D_HP_ANISO_H = c_H2D_HP_ANISO_H
    H2D_HP_ANISO_P = c_H2D_HP_ANISO_P
    H2D_HP_ANISO = c_H2D_HP_ANISO

cdef class SelOption:
    """ Options of the selector.

        Possible values are attributes of the class:

    - PREFER_SYMMETRIC_MESH: Prefer symmetric mesh when selection of the best candidate is done. If two or more candiates has the same score, they are skipped. This option is set by default.
    - APPLY_CONV_EXP_DOF: Use $d^c - d_0^c$, where $c$ is the convergence exponent, instead of $(d - d_0)^c$ to evaluate the score. This option is not set by default.
    """
    PREFER_SYMMETRIC_MESH = c_H2D_PREFER_SYMMETRIC_MESH
    APPLY_CONV_EXP_DOF = c_H2D_APPLY_CONV_EXP_DOF

cdef class ProjBasedSelector:
    """ A base class for projection based selector.

        Suggests a refinement based on en error of proposed refinements (candidates).
    An error of a candidate is a combination of errors of
    elements of a candidate. Each element of a candidate is calculated
    separatelly.

        Do not use this class directly,
    use derived classes H1ProjBasedSelector and L2ProjBasedSelector instead.
    """
    def __cinit__(self):
        pass 

    def __dealloc__(self):
        delete(self.thisptr)

    def set_error_weights(self, double weight_h, double weight_p, double weight_aniso):
        """ Sets error weights.
        
        An error weight is a multiplicative coefficient that modifies an error of candidate.
        Default wieghts are:

        - H-candidate weight: 2.0
        - P-candidate weight: 1.0
        - ANISO-candidate weight: 1.4142

        Suggested Use:

        - settings the default weights: ```selector.set_error_weights(2.0, 1.0, 1.4142)```
        """
        self.thisptr.set_error_weights(weight_h, weight_p, weight_aniso)

    def set_option(self, int option, int enable):
        """ Enables or disables an option.

        Suggested Use:

        - enabling an option ``PREFER_SYMMETRIC_MESH``: ``selector.set_option(SelOption.PREFER_SYMMETRIC_MESH, true)``
        - disabling an option ``APPLY_CONV_EXP_DOF``: ``selector.set_optin(SelOption.APPLY_CONV_EXP_DOF, false)``
    """
        self.thisptr.set_option(<c_SelOption>option, enable)

cdef class H1ProjBasedSelector(ProjBasedSelector):
    """ A projection based selector for H1 space by the method adapt() of the class Adapt
        The selector should be instantiated outside the adaptivity loop.

    Suggested Use:

    ```selector = H1ProjBasedSelector(cand_list, conv_exp, max_order, h1_shapeset)```

    where:

    - ```cand_list``` defines possible candidates, e.g., CandList.HP_ANISO.
    - ```conv_exp``` specifies a convergence exponent which incluences selecting of a candidate. For more detials, refer to User Documentation.
    - ```max_order``` defined a maximum order which should be used by the refinement. Currently, the maximum order is 9. Supply -1 to use the maximum supported order.
    - ```h1_shapeset``` expects an instance of the class H1Shapeset.
    """
    def __cinit__(self, int cand_list, double conv_exp, int max_order, H1Shapeset shapeset = None):
        self.thisptr = <c_ProjBasedSelector*>(new_H1ProjBasedSelector(<c_CandList>cand_list, conv_exp, max_order, NULL))

cdef class L2ProjBasedSelector(ProjBasedSelector):
    """ A projection based selector for L2 space used by the method adapt() of the class Adapt.
        The selector should be instantiated outside the adaptivity loop.

    Suggested Use:

        ```selector = L2ProjBasedSelector(cand_list, conv_exp, max_order, l2_shapeset)```

    where:

    - ```cand_list``` defines possible candidates, e.g., CandList.HP_ANISO.
    - ```conv_exp``` specifies a convergence exponent which incluences selecting of a candidate. For more detials, refer to User Documentation.
    - ```max_order``` defined a maximum order which should be used by the refinement. Currently, the maximum order is 9. Supply -1 to use the maximum supported order.
    - ```l2_shapeset``` expects an instance of the class L2Shapeset.
    """
    def __cinit__(self, int cand_list, double conv_exp, int max_order, L2Shapeset shapeset):
        self.thisptr = <c_ProjBasedSelector*>(new_L2ProjBasedSelector(<c_CandList>cand_list, conv_exp, max_order, shapeset.thisptr))

cdef class Adapt:
    """ A base class for adaptivity. Adaptivity provides framework for modyfying elements in order to decrease errors of the solution.
        Do not use this class directly, use derived classes H1Adapt and L2Adapt instead."""
    def __cinit__(self):
        pass

    def __dealloc__(self):
        delete(self.thisptr)

    def set_solutions(self, sln_list, rsln_list):
        """ Sets coarse solutions and reference solutions.
        This method has to be called before the method calc_error().

        Suggested Use:

        - single component: ```hp.set_solutions([sln1], [rsln1])```
        - multiple components: ```hp.set_solutions([sln1, sln2], [rsl1, rsln2])```
    """
        cdef c_SolutionTuple slns
        cdef c_SolutionTuple rslns
        for i in range(len(sln_list)):
            slns.push_back(<c_Solution*>((<Solution>(sln_list[i])).thisptr))
        for i in range(len(rsln_list)):
            rslns.push_back(<c_Solution*>((<Solution>(rsln_list[i])).thisptr))
        self.thisptr.set_solutions(slns, rslns)

    def calc_error(self, int error_flags = H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS):
        """ Calculates errors of elements and returns a total error.
        This method has to be called before the method adapt().

        Suggested Use: ```err_percent = hp.calc_error() * 100```
    """
        return self.thisptr.calc_error(error_flags)

    def adapt(self, ProjBasedSelector selector, double thr, int strat = 0, int regularize = -1, double to_be_processed = 0.0):
        """ Does adaptive refinement based on errors of elements.
        Refinements are selected using a supplied selector.

        Returns ```True``` if no element was refined.

        Suggested Use:
        
            ```hp.adapt(selector, thr, strat, regularize, same_order, to_be_processed)```

        where

        - ```selector``` is a selector, see the class ProjBasedSelector.
        - ```thr``` is a threshold used by adaptivity strategy.
        - ```strat``` is an adaptivity strategy. Possible values are 0, 1, 2, and 3.
        - ```regularize``` specifies mesh regularization. Default value is `-1`.
        - ```same_orders``` specifies whether all element should have the same order. Default is `False`.
        - ```to_be_procesed``` specifies an error to process. Used by strategy 3. Default is `0.0`.
        """
        return self.thisptr.adapt(selector.thisptr, thr, strat, regularize, to_be_processed)

cdef class H1Adapt(Adapt):
    """ Adaptivity class for H1 space.
        For details on class members, see the parent class Adapt.

        Suggested Use:

    - single component: ```hp = H1Adapt([h1_space1])```
    - multiple components: ```hp = H1Adapt([h1_space1, h1_space2])```
    """
    """
    def __cinit__(self, space_list):
        cdef c_H1SpaceTuple spaces
        for i in range(len(space_list)):
            spaces.push_back((<H1Space>(space_list[i])).thisptr)
        self.thisptr = <c_Adapt*>(new_H1Adapt(spaces))
    """
    
    def __cinit__(self, LinSystem ls):
        self.thisptr = <c_Adapt*>(new_H1Adapt(ls.thisptr))
            
cdef class L2Adapt(Adapt):
    """ Adaptivity class for L2 space.
        For details on class members, see the parent class Adapt.

        Suggested Use:

    - single component: ```hp = L2Adapt([l2_space1])```
    - multiple components: ```hp = L2Adapt([l2_space1, l2_space2])```
    """
    """
    def __cinit__(self, space_list):
        cdef c_L2SpaceTuple spaces
        for i in range(len(space_list)):
            spaces.push_back((<L2Space>(space_list[i])).thisptr)
        self.thisptr = <c_Adapt*>(new_L2Adapt(spaces))
    """
    def __cinit__(self, LinSystem ls):
        self.thisptr = <c_Adapt*>(new_L2Adapt(ls.thisptr))    

cdef class Linearizer:
    """
    Linearizes the solution.

    It returns the triangles and vertices and you can then use it to visualize
    the solution.

    Example::

        In [40]: l = Linearizer(sln)

        In [44]: l.process_solution(sln)
        Linearizer: 5519 verts, 10847 tris in 0.05 sec

        In [45]: l.get_vertices()
        Out[45]:
        array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
               [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
               [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
               ...,
               [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
               [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
               [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])

        In [46]: l.get_triangles()
        Out[46]:
        array([[   3, 5448,   29],
               [  27, 5445,   28],
               [  29,   28,   26],
               ...,
               [5499, 5498, 5479],
               [5510, 5493, 5491],
               [5513, 5508, 5491]], dtype=int32)

    """

    def __cinit__(self):
        self.thisptr = new_Linearizer()

    def __dealloc__(self):
        delete(self.thisptr)

    def process_solution(self, MeshFunction sln):
        self.thisptr.process_solution(<c_MeshFunction *>sln.thisptr)

    def get_vertices(self):
        """
        Returns the list of vertices.

        It's a list of triples, where each triple contains (x, y, val), where
        x, y are the "x" and "y" coordinates of the vertex and "val" is the
        value of the solution at the vertex.

        Example::

            In [45]: l.get_vertices()
            Out[45]:
            array([[  0.00000000e+00,  -1.00000000e+00,  -2.22396971e-17],
                   [  1.00000000e+00,  -1.00000000e+00,  -1.64798730e-17],
                   [ -1.00000000e+00,   0.00000000e+00,   8.09899023e-17],
                   ...,
                   [  1.48437500e-01,  -1.56250000e-02,   1.62359362e-01],
                   [  1.32812500e-01,   0.00000000e+00,   1.56012622e-01],
                   [  1.32812500e-01,  -1.56250000e-02,   1.50562411e-01]])

        """
        cdef double3 *vert = self.thisptr.get_vertices()
        cdef int nvert = self.thisptr.get_num_vertices()
        cdef ndarray vec = c2numpy_double(<double *>vert, 3*nvert)
        return vec.reshape((nvert, 3))

    def get_num_vertices(self):
        return self.thisptr.get_num_vertices()

    def get_triangles(self):
        """
        Returns a list of triangles.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

            In [46]: l.get_triangles()
            Out[46]:
            array([[   3, 5448,   29],
                   [  27, 5445,   28],
                   [  29,   28,   26],
                   ...,
                   [5499, 5498, 5479],
                   [5510, 5493, 5491],
                   [5513, 5508, 5491]], dtype=int32)

        """
        cdef int3 *tri = self.thisptr.get_triangles()
        cdef int ntri = self.thisptr.get_num_triangles()
        cdef ndarray vec = c2numpy_int(<int *>tri, 3*ntri)
        return vec.reshape((ntri, 3))

    def get_num_triangles(self):
        return self.thisptr.get_num_triangles()

    def get_edges(self):
        """
        Returns a list of edges.

        The list contains triples of vertices IDs. Use get_vertices() to obtain
        vertices coordinates.

        Example::

            In [47]: l.get_edges()
            Out[47]:
            array([[   3,   27,    0],
                   [  27,   24,    0],
                   [  24,   30,    0],
                   ...,
                   [5339, 5070,    4],
                   [5070, 5077,    4],
                   [5077,   11,    4]], dtype=int32)

        """
        cdef int3 *edges = self.thisptr.get_edges()
        cdef int nedges = self.thisptr.get_num_edges()
        cdef ndarray vec = c2numpy_int(<int *>edges, 3*nedges)
        return vec.reshape((nedges, 3))

    def get_num_edges(self):
        return self.thisptr.get_num_edges()

    def get_min_value(self):
        return self.thisptr.get_min_value()

    def get_max_value(self):
        return self.thisptr.get_max_value()

    def save_data(self, char* filename):
        self.thisptr.save_data(filename)

    def load_data(self, char* filename):
        self.thisptr.load_data(filename)

cdef class Vectorizer(Linearizer):

    def __cinit__(self):
        self.thisptr = <c_Linearizer *>new_Vectorizer()

    # this segfaults:
    #def __dealloc__(self):
    #    delete(<c_Vectorizer *>(self.thisptr))

    def process_solution(self, MeshFunction xsln, MeshFunction ysln,
            int xitem=c_H2D_FN_VAL_0, int yitem=c_H2D_FN_VAL_0, double eps=c_H2D_EPS_LOW):
        (<c_Vectorizer *>(self.thisptr)).process_solution(
                <c_MeshFunction *>xsln.thisptr, xitem,
                <c_MeshFunction *>ysln.thisptr, yitem,
                eps)

    def get_vertices(self):
        cdef c_Vectorizer *_self = <c_Vectorizer *>(self.thisptr)
        cdef double4 *vert = _self.get_vertices()
        cdef int nvert = self.thisptr.get_num_vertices()
        cdef ndarray vec = c2numpy_double(<double *>vert, 4*nvert)
        return vec.reshape((nvert, 4))

    def get_dashes(self):
        cdef c_Vectorizer *_self = <c_Vectorizer *>(self.thisptr)
        cdef int2 *dashes = _self.get_dashes()
        cdef int ndashes = _self.get_num_dashes()
        cdef ndarray vec = c2numpy_int(<int *>dashes, 2*ndashes)
        return vec.reshape((ndashes, 2))

cdef class View:
    def wait(self):
        View_wait()

cdef class ScalarView(View):
    cdef c_ScalarView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
                int height = 800):
        self.thisptr = new_ScalarView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, MeshFunction sln, eps=None):
        cdef Function s = <Function>sln
        if eps is None:
            self.thisptr.show(<c_MeshFunction *>s.thisptr)
        else:
            self.thisptr.show(<c_MeshFunction *>s.thisptr, <double>eps)

    def show_mesh(self, show=True):
        self.thisptr.show_mesh(show)

    def show_scale(self, show=True):
        self.thisptr.show_scale(show)

    def set_min_max_range(self, double min, double max):
        self.thisptr.set_min_max_range(min, max)

    def set_title(self, char *title):
        ((self.thisptr)).set_title(title)
        #(<c_View *>(self.thisptr)).set_title(title)

    def wait(self):
        self.thisptr.wait()

cdef class BaseView(View):
    cdef c_BaseView *thisptr

    def __cinit__(self, char *title="BaseView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_BaseView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class MeshView(View):
    cdef c_MeshView *thisptr

    def __cinit__(self, char *title="MeshView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_MeshView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, Mesh s):
        self.thisptr.show(s.thisptr)

cdef class OrderView(View):
    cdef c_OrderView *thisptr

    def __cinit__(self, char *title="OrderView", int x=-1, int y=-1,
            int width=1000, int height = 800):
        self.thisptr = new_OrderView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, H1Space s):
        self.thisptr.show(s.thisptr)

cdef class VectorView(View):
    cdef c_VectorView *thisptr

    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
            int height = 800):
        self.thisptr = new_VectorView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    def show(self, MeshFunction s1, MeshFunction s2, eps = 0.008):
        (<c_VectorView *>(self.thisptr)).show(<c_MeshFunction *>s1.thisptr, <c_MeshFunction *>s2.thisptr, eps)

    def show_scale(self, show=True):
        self.thisptr.show_scale(show)

    def set_min_max_range(self, double min, double max):
        self.thisptr.set_min_max_range(min, max)

#cdef class MatrixView(View):
#    cdef c_MatrixView *thisptr
#
#    def __cinit__(self, char *title, int x=-1, int y=-1, int width=1000,
#            int height = 800):
#        self.thisptr = new_MatrixView(title, x, y, width, height)

    #def __dealloc__(self):
    #    delete(self.thisptr)

    #def show(self, DiscreteProblem ep):
    #    self.thisptr.show(ep.thisptr)

def init_hermes2d_wrappers():
    init_global_empty_tuple()

def set_verbose(verbose):
    """
    Sets the global verbose_mode variable.

    That variable controls how verbose hermes2d is.

    Returns the old status.
    """
    global c_verbose_mode
    global c_info_mode
    old_flag = c_verbose_mode
    c_verbose_mode = verbose
    c_info_mode = verbose
    return old_flag

def set_warn_integration(warn_integration):
    """
    Sets the global warn_integration variable.

    That variable controls if warn_order() in the C++ hermes2d should emit the
    warning about integration rules.

    Returns the old status.
    """
    global c_warn_integration
    old_flag = c_warn_integration
    c_warn_integration = warn_integration
    return old_flag

#def glut_main_loop():
#    """
#    This function waits for all Views to finish.
#
#    E.g. until the user closes them. If you don't call this function, all
#    windows will be closed when the program ends, so you need to call it if you
#    want to give the user a chance to manipulate with the plots at the end of
#    your script (program).
#    """
#    finish_glut_main_loop(False)

init_hermes2d_wrappers()
