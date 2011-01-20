import math

class ParseError(Exception):
    pass

def convert2tuple(s):
    """
    Converts any iterable to a tuple recursively.

    Insert your iterable into parameter "s".

    Example:

    >>> convert2tuple([[0,0],[0,1],[1,1],[1,0],[0.25,0.25],[0.25,0.75],[0.75,0.75],[0.75,0.25]]) 
    ((0, 0), (0, 1), (1, 1), (1, 0), (0.25, 0.25), (0.25, 0.75), (0.75,
    0.75), (0.75, 0.25)) 
    >>> convert2tuple([[0,0],[0,1],[1,1],[1,0],[0.25,0.25],[0.25,0.75],[0.75,0.5]]) 
    ((0, 0), (0, 1), (1, 1), (1, 0), (0.25, 0.25), (0.25, 0.75), (0.75,
    0.5)) 

    """
    if hasattr(s, "__iter__"):
        return tuple([convert2tuple(y) for y in s])
    else:
        return s

def read_hermes_format(filename):
    """
    Reads a mesh from a file in a hermes format.

    Returns nodes, elements, boundary, nurbs or raises a ParseError if the
    syntax is invalid.
    """
    m = open(filename).read()
    return read_hermes_format_str(m)

def read_hermes_format_str(m):
    """
    Reads a mesh from a string in a hermes format.

    Returns nodes, elements, boundary, nurbs or raises a ParseError if the
    syntax is invalid.
    """
    m = m.strip()
    m = m.replace("\r\n", "\n")
    m = m.replace("=\n", "= \\\n")
    m = m.replace("{", "[")
    m = m.replace("}", "]")
    m = m.replace("^", "**")
    # Make sure 1/2 produces 0.5:
    m = "from __future__ import division\n\n" + m
    namespace = {}
    try:
        exec m in math.__dict__, namespace
    except (SyntaxError, NameError), e:
        print "Error while parsing the mesh file. More information below."
        raise ParseError(str(e))
    nodes = namespace.pop("vertices", None)
    elements = namespace.pop("elements", None)
    boundary = namespace.pop("boundaries", None)
    nurbs = namespace.pop("curves", None)
    refinements = namespace.pop("refinements", None)
    if nodes is None or elements is None or boundary is None:
        raise ParseError("Either nodes, elements or boundary is missing")
    return convert2tuple(nodes), convert2tuple(elements), \
            convert2tuple(boundary), convert2tuple(nurbs), \
            convert2tuple(refinements)
