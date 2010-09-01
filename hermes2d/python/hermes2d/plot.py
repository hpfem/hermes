from hermes2d import Linearizer, Solution

def sln2png(sln, filename):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln)
    from enthought.mayavi.mlab import savefig
    savefig(filename)

def plot_sln_mpl(sln, method="default", just_mesh=False, axes=None):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    method = "default" ... creates a plot using triangles (the triangles are
                not interpolated, so sometimes one can see small defects)
    method = "contour" ... takes the vertices from linearizer and interpolates
                them using contour and contourf (it doesn't take into account
                the triangulation, so one can see the defects from the convex
                hull approximation)

    just_mesh ... only shows the mesh, but not the solution
    """
    lin = Linearizer()
    lin.process_solution(sln)
    v = lin.get_vertices()
    if method=="contour":
        from numpy import min, max, linspace
        from matplotlib.mlab import griddata
        import matplotlib.pyplot as plt
        x = v[:, 0]
        y = v[:, 1]
        z = v[:, 2]
        # define grid.
        xi = linspace(min(x), max(x), 100)
        yi = linspace(min(y), max(y), 100)
        # grid the data.
        zi = griddata(x, y, z, xi, yi)
        # contour the gridded data, plotting dots at the nonuniform data points.
        CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
        CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('Solution')
    elif method == "default":
        from numpy import array
        import matplotlib.collections as collections
        #import matplotlib.pyplot as plt
        if axes is None:
            from pylab import gca
            axes = gca()
        verts = []
        vals = []
        for t in lin.get_triangles():
            triangle = tuple([tuple(v[n][:2]) for n in t])
            val = sum([v[n][2] for n in t])
            vals.append(val/3.)
            verts.append(triangle)
        verts = array(verts)
        vals = array(vals)
        if just_mesh:
            lw = 1
        else:
            lw = 0
        col = collections.PolyCollection(verts, linewidths=lw, antialiaseds=0)
        col.set_array(vals)
        #col.set_cmap(plt.cm.jet)
        ax = axes
        ax.add_collection(col)
        ax.set_xlim(verts[:, :, 0].min(), verts[:, :, 0].max())
        ax.set_ylim(verts[:, :, 1].min(), verts[:, :, 1].max())
        ax.set_aspect("equal")
        #plt.colorbar()
        #plt.title('Solution')
    else:
        raise ValueError("Unknown method (%s)" % method)

def plot_sln_mayavi(sln, notebook=False):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    Currently only a very simple version is implemented, that takes the
    vertices from linearizer and interpolates them. More sophisticated version
    should take the triangles.
    """
    lin = Linearizer()
    lin.process_solution(sln)
    vert = lin.get_vertices()
    triangles = lin.get_triangles()
    from numpy import zeros
    from enthought.mayavi import mlab
    x = vert[:, 0]
    y = vert[:, 1]
    z = zeros(len(y))
    t = vert[:, 2]
    if notebook:
        # the off screen rendering properly works only with VTK-5.2 or above:
        mlab.options.offscreen = True
    mlab.clf()
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    mlab.view(0, 0)

    # Below is a code that does exactly what the "View along the +Z axis"
    # button does:
    #scene = mlab.get_engine().current_scene.scene
    #scene.camera.focal_point = [0, 0, 0]
    #scene.camera.position = [0, 0, 1]
    #scene.camera.view_up = [0, 1, 0]
    #scene.renderer.reset_camera()
    #scene.render()
    # the above looks ok, but there is still quite a large margin, so we prefer
    # to just call .view(0, 0), which seems to be working fine.
    return s

def plot_hermes_mesh_mpl(mesh, space=None, edges_only=False):
    if space is None:
        polynomial_orders = None
    else:
        polynomial_orders = mesh.get_elements_order(space)
    return plot_mesh_mpl(mesh.nodes_dict, mesh.elements,
                    mesh.get_polygonal_boundary(),
                    polynomial_orders=polynomial_orders,
                    edges_only=edges_only)

def plot_mesh_mpl(nodes, elements, polygons=None,
        polynomial_orders=None, edges_only=False):
    from matplotlib import pyplot
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    from matplotlib.patches import Rectangle
    colors_old = {
            1: '#000684',
            2: '#3250fc',
            3: '#36c4ee',
            4: '#04eabc',
            5: '#62ff2a',
            6: '#fdff07',
            7: '#ffa044',
            8: '#ff1111',
            9: '#b02c2c',
            10: '#820f97',
            }
    colors = {
            0: '#7f7f7f',
            1: '#7f2aff',
            2: '#2a2aff',
            3: '#2a7fff',
            4: '#00d4aa',
            5: '#00aa44',
            6: '#abc837',
            7: '#ffd42a',
            8: '#c87137',
            9: '#c83737',
            10: '#ff0000',
            }
    fig = pyplot.figure()
    sp = fig.add_subplot(111)
    for el_id in polygons:
        x = list(polygons[el_id][:, 0])
        y = list(polygons[el_id][:, 1])
        x.append(x[0])
        y.append(y[0])
        vertices = zip(x, y)
        codes = [Path.MOVETO] + [Path.LINETO]*(len(vertices)-2) + \
                    [Path.CLOSEPOLY]
        p = Path(vertices, codes)
        if edges_only:
            color = "white"
            linewidth = 2
        else:
            if polynomial_orders is None:
                color = colors[0]
            else:
                color = colors[polynomial_orders[el_id]]
            linewidth = 1
        patch = PathPatch(p, facecolor=color, lw=linewidth,
                edgecolor='#000000')
        sp.add_patch(patch)
    show_legend = polynomial_orders is not None

    if show_legend:
        # Create legend
        def split_nodes():
            x = []
            y = []

            if isinstance(nodes, dict):
                _nodes = nodes.items()
            else:
                _nodes = enumerate(nodes)
            for k, pnt in _nodes:
                x.append(pnt[0])
                y.append(pnt[1])

            return (x, y)

        def get_max(what='x'):
            x, y = split_nodes()

            if what == 'x':
                return max(x)
            else:
                return max(y)

        def get_min(what='x'):
            x, y = split_nodes()

            if what == 'x':
                return min(x)
            else:
                return min(y)

        maxX = get_max('x')
        maxY = get_max('y')

        minX = get_min('x')
        minY = get_min('y')

        dy = (maxY - minY) / 20
        dx = (maxX - minX) / 20

        y = minY + dy
        x = maxX + dx

        ord = polynomial_orders.items()
        order_list = []
        for k,v in ord:
            order_list.append(v)
        m = max(order_list)

        for k,c in colors.items():
            if k <= m :
                p = Rectangle(xy=(x,y), width=dx, height=dy, fill=True, facecolor=c)
                sp.add_patch(p)
                sp.text(x + dx + (dx/2), y + (dy/4), str(k))
                y += dy
            else:
                break

        sp.text(x, y + (dy/2), str('Orders'))
    sp.set_title("Mesh")
    sp.set_aspect("equal")
    sp.autoscale_view()
    return sp.figure

class ScalarView(object):

    def __init__(self, name="Solution", x=0, y=0, w=50, h=50):
        self._name = name
        self._lib = None
        self._notebook = False
        self._glut_view = None

    def show_scale(self, *args):
        pass

    def show_mesh(self, *args):
        pass

    def wait(self):
        if self._lib == "mpl" and self._notebook == False:
            import pylab
            pylab.show()

    def show(self, sln, show=True, lib="mayavi", notebook=None,
            filename="scalar.png", **options):
        """
        Shows the solution.

        show ... should it actually plot the window? Set to False in tests.
        lib .... which library to use for the plotting? either "mpl" or "mayavi"
        notebook ... are we running inside Sage notebook? If True, just save
                the image to a.png
        filename ... the name of the filename if we are saving the image (e.g.
                notebook == False)

        Example:

        >>> 1 + 1
        2
        >>> 1 + 2
        3
        """
        if notebook is None:
            try:
                from sagenb.misc.support import EMBEDDED_MODE
            except ImportError:
                EMBEDDED_MODE = False
            notebook = EMBEDDED_MODE

        self._lib = lib
        self._notebook = notebook
        if lib == "glut":
            if self._glut_view is None:
                from _hermes2d import ScalarView
                self._glut_view = ScalarView(self._name)
            self._glut_view.show(sln)
        elif lib == "mpl":
            plot_sln_mpl(sln, **options)
            import pylab
            if show:
                if notebook:
                    pylab.savefig(filename)
                else:
                    pylab.ion()
                    pylab.draw()
                    pylab.ioff()
        elif lib == "mayavi":
            plot_sln_mayavi(sln, notebook=notebook)
            from enthought.mayavi import mlab
            if show:
                engine = mlab.get_engine()
                image = engine.current_scene
                image.scene.background = (1.0, 1.0, 1.0)
                image.scene.foreground = (0.0, 0.0, 0.0)
                #mlab.colorbar(orientation="vertical")
                if notebook:
                    mlab.savefig(filename)
                else:
                    mlab.show()
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)

class MeshView(object):

    def __init__(self, name="Solution", x=0, y=0, w=500, h=500):
        self._name = name
        self._x = x
        self._y = y
        self._w = w
        self._h = h

    def wait(self):
        pass

    def show(self, mesh, show=True, lib="mpl", notebook=None, space=None,
            filename="mesh.png", **options):
        if notebook is None:
            try:
                from sagenb.misc.support import EMBEDDED_MODE
            except ImportError:
                EMBEDDED_MODE = False
            notebook = EMBEDDED_MODE
        if lib == "glut":
            from _hermes2d import MeshView
            m = MeshView(self._name, self._x, self._y, self._w, self._h)
            m.show(mesh)
            m.wait()
        elif lib == "mpl":
            p = plot_hermes_mesh_mpl(mesh, space=space, **options)
            if show:
                if notebook:
                    p.savefig(filename)
                else:
                    p.show()
                    import pylab
                    pylab.show()
            return p
        else:
            raise NotImplementedError("Unknown library '%s'" % lib)
