===================
Welcome to Hermes3D
===================

Hermes3D is a C++ library for rapid prototyping of adaptive hp-FEM solvers for
3D problems. The internal structure, API, and usage of Hermes3D is very similar
to Hermes2D. Please get familiar with Hermes2D before working with Hermes3D. 
The User Documentation and Tutorial follow the same philosophy.


License
=======

Hermes2D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public
License along with Hermes2D. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA  02110-1301, USA.


Tests (specific details for Hermes3D)
=====================================

To enable tests, say 'set(WITH_TESTS YES)' in your CMake.vars.
To run the tests, type::

    $ make test

To run quick tests, type::

    $ make test-quick

Note: To run developer tests, say 'set(DEV_TESTS YES)' in CMake.vars. This is
needed to run only if the lowest internals are changed. Developers test suite
includes hundreds of tests for hanging nodes. These will run for several hours
in case of H1 space, the Hcurl ones take days. These test do not have to be run
every time.
