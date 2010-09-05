"""
This module contains testing classes for the leaks test.
"""

class A(object):

    def __init__(self, c):
        print "A.__init__() was called"
        self._c = c

    def __del__(self):
        print "A.__del__() was called"
        self._c.event("__del__ was called")
