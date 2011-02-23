Multimesh hp-FEM
----------------

In multiphysics PDE problems (and other PDE systems) it often happens that one
physical field (solution component) has a singularity or a boundary layer 
where other fields are smooth. If one approximates all of them on the 
same mesh, then some of them will be refined where this is absolutely not needed.
This can be extremely wasteful, as we will see in the tutorial example 02-system-adapt. 
But first let us introduce the multimesh discretization method that we developed 
to circumvent this problem.

Multimesh assembling
~~~~~~~~~~~~~~~~~~~~

Hermes makes it possible to approximate each physical field (or solution
component) on an individual mesh. These meshes are not completely independent
of each other -- they have a common coarse mesh that we call *master mesh*.
The master mesh is there for algorithmic purposes only, it may not 
even be used for discretization purposes: Every mesh in the system 
is obtained from it via an arbitrary sequence of elementary refinements.
This is illustrated in the following figure, where (A) is the master mesh,
(B) - (D) three different meshes (say, for a coupled problem with three
equations), and (E) is the virtual *union mesh* that is used for assembling.

.. image:: multimesh-example/multimesh.png
   :align: center
   :width: 750
   :alt: Multimesh

The union mesh is not constructed physically in the computer memory -- 
merely it serves as a hint to correctly transform integration points
while integrating over sub-elements of the elements of the existing meshes. 
The following figure shows the integration over an element $Q_k$ of the 
virtual union mesh, and what are the appropriate subelements of the 
existing elements where this integration is performed:

.. image:: multimesh-example/multimesh2.png
   :align: center
   :width: 600
   :alt: Multimesh

As a result, the multimesh discretization of the PDE system is *monolithic*
in the sense that *no physics is lost* -- all integrals in the 
discrete weak formulations are evaluated exactly up to the error in the 
numerical quadrature. In particular, we do not perform operator splitting 
or commit errors while transferring solution data between different meshes.
The multimesh assembling in Hermes works with all meshes at the same time, 
there is no such thing as interpolating or projecting functions between 
different meshes. More details about this method can be found in the 
corresponding `research articles <http://hpfem.org/hermes/doc/src/citing-hermes.html>`_. 

Simultaneous adaptivity on multiple meshes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As explained above, each physical field (solution component0 has its own mesh 
and these meshes can be very different. The adaptivity algorithm puts all elements 
of all these meshes one single array, and calculates an error function as a difference
between the coarse and reference mesh approximations. Elements are then sorted according 
to their norm of the error function, and then the ones with the largest errors are refined
as in the standard hp-FEM. Note that if some physical field is already resolved well, 
this algorithm leaves its mesh alone and will only perform refinements in the other 
meshes where they are needed. 

The function calc_error()
~~~~~~~~~~~~~~~~~~~~~~~~~

Errors for the adaptvity algorithm as well as the total error are calculated in the 
function calc_error().
This function can be called with four combinations of the absolute/relative flags 
for the total error and the element errors::

    hp.calc_error(H2D_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL)
    hp.calc_error(H2D_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS)
    hp.calc_error(H2D_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_REL)
    hp.calc_error(H2D_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS)

The exact meaning of these flags is as follows:

- ```H2D_TOTAL_ERROR_REL```: Total error is the norm of the error function (which is the difference between the reference and coarse mesh solutions) divided by the norm of the reference solution.
- ```H2D_TOTAL_ERROR_ABS```: Total error is the norm of the error function.
- ```H2D_TOTAL_ERROR_REL```: Element error is the norm of the error function on that element divided by the norm of the corresponding solution component. 
- ```H2D_TOTAL_ERROR_ABS```: Element error is the norm of the error function on that element.

