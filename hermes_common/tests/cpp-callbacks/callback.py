"""
This module contains testing classes for the cpp-callbacks test.
"""

class A(object):

    def __init__(self):
        print "A.__init__() was called"

    def __del__(self):
        print "A.__del__() was called"
