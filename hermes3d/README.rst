Specific Hermes3D instructions
==============================

Read the main README file first (one directory up). This file contains specific
instructions for Hermes3D, that you can't find in the main README file.

Hermes3D is a C++ library for rapid prototyping of adaptive hp-FEM solvers for
3D problems. The usage of Hermes3D is very similar to Hermes2D, and thus it is
recommend to get familiar with Hermes2D first. Also the User Documentation
of Hermes3D concentrates on things that are different from Hermes2D.

Tests
-----

To enable tests, say 'set(WITH_TESTS YES)' in your CMake.vars.
To run the tests, type::

    $ make test

To run quick tests, type::

    $ make test-quick

Note: To run developer tests, say 'set(DEV_TESTS YES)' in CMake.vars. This is
needed to run only if the lowest internals are changed. Developers test suite
includes hundreds of tests for hanging nodes. These will run for several hours
in case of H1 space, Hcurl ones take days. These test do not have to be run
every time.
