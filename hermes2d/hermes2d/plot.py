from hermes2d import Linearizer, Solution

def sln2png(sln, filename, offscreen=False):
    """
    Creates a nice png image of the Solution sln.
    """
    plot_sln_mayavi(sln, offscreen=offscreen)
    from enthought.mayavi.mlab import savefig
    savefig(filename)


def plot_sln_mayavi(sln, offscreen=False, show_scale=True):
    """
    Plots the Solution() instance sln using Linearizer() and matplotlib.

    It takes the vertices from linearizer and interpolates them.
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
    if offscreen:
        # the off screen rendering properly works only with VTK-5.2 or above:
        mlab.options.offscreen = True
    mlab.clf()
    mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    mlab.view(0, 0)
    mlab.view(distance=4)
    mlab.view(focalpoint=(.35,0,0))
    mlab.colorbar(title="Solution", orientation="vertical")

    #mlab.move(right=-1.0, up=-10.0)
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
    return mlab
