Nernst-Planck Equation System
-----------------------------

**Git reference:** Example `nernst-planck-timedep-adapt <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/examples/nernst-planck/nernst-planck-timedep-adapt>`_.

**Equation reference:** The first version of the following derivation was published in:
*IPMC: recent progress in modeling, manufacturing, and new applications 
D. Pugal, S. J. Kim, K. J. Kim, and K. K. Leang 
Proc. SPIE 7642, (2010)*.
The following Bibtex entry can be used for the reference:

::

	@conference{pugal:76420U,
		author = {D. Pugal and S. J. Kim and K. J. Kim and K. K. Leang},
		editor = {Yoseph Bar-Cohen},
		title = {IPMC: recent progress in modeling, manufacturing, and new applications},
		publisher = {SPIE},
		year = {2010},
		journal = {Electroactive Polymer Actuators and Devices (EAPAD) 2010},
		volume = {7642},
		number = {1},
		numpages = {10},
		pages = {76420U},
		location = {San Diego, CA, USA},
		url = {http://link.aip.org/link/?PSI/7642/76420U/1},
		doi = {10.1117/12.848281}
	}

The example is concerned with the finite element solution 
of the Poisson and Nernst-Planck equation system. The Nernst-Planck
equation is often used to describe the diffusion, convection,
and migration of charged particles:

.. math::
	:label: nernstplanck

		\frac {\partial C} {\partial t} + \nabla \cdot 
		(-D\nabla C - z \mu F C \nabla \phi) = 
		- \vec {u} \cdot \nabla C.

The second term on the left side is diffusion and the third term is
the migration that is directly related to the the local voltage
(often externally applied) $\phi$. The term on the right side is
convection. This is not considered in the current example. The variable
$C$ is the concentration of the particles at any point of a domain
and this is the unknown of the equation.

One application for the equation is to calculate charge configuration
in ionic polymer transducers. Ionic polymer-metal composite is
for instance an electromechanical actuator which is basically a thin
polymer sheet that is coated with precious metal electrodes on both
sides. The polymer contains fixed anions and mobile cations such
as $H^{+}$, $Na^{+}$ along with some kind of solvent, most often water.

When an voltage $V$ is applied to the electrodes, the mobile cations
start to migrate whereas immobile anions remain attached to the polymer
backbone. This creates spatial charges, especially near the electrodes.
One way to describe this system is to solve Nernst-Planck equation
for mobile cations and use Poisson equation to describe the electric
field formation inside the polymer. The poisson equation is

.. math::
	:label: poisson

		\nabla \cdot \vec{E} = \frac{F \cdot \rho}{\varepsilon},

where $E$ could be written as $\nabla \phi = - \vec{E}$ and $\rho$ is
charge density, $F$ is the Faraday constant and $\varepsilon$ is dielectric
permittivity. The term $\rho$ could be written as:

.. math::
	:label: rho
	
		\rho = C - C_{anion},
		
where $C_{anion}$ is a constant and equals anion concentration. Apparently
for IPMC, the initial spatial concentration of anions and cations are equal.
The inital configuration is shown:

.. image:: example-np/IPMC.png
	:align: center
	:width: 377
	:height: 173
	:alt: Initial configuration of IPMC.

The purploe dots are mobile cations. When a voltage is applied, the anions
drift:

.. image:: example-np/IPMC_bent.png
	:align: center
	:width: 385
	:height: 290
	:alt: Bent IPMC

Images reference: 
*IPMC: recent progress in modeling, manufacturing, and new applications 
D. Pugal, S. J. Kim, K. J. Kim, and K. K. Leang 
Proc. SPIE 7642, (2010)*
This eventually results in actuation (mostly bending) of the material (not considered in this section).

To solve equations :eq:`nernstplanck` and :eq:`poisson` boundary conditions must be specified as well.
When solving in 2D, just a cross section is considered. The boundaries are
shown in: 

.. image:: example-np/IPMC_schematic.png
	:align: center
	:width: 409 
	:height: 140
	:alt: IPMC boundaries

For Nernst-Planck equation :eq:`nernstplanck`, all the boundaries have the same, insulation
boundary conditions:

.. math::
	:label: nernstboundary

	-D \frac{\partial C}{\partial n} - z \mu F C \frac{\partial \phi} {\partial n} = 0

For Poisson equation:

 #. (positive voltage): Dirichlet boundary $\phi = 1V$. For some cases it might be necessary to use electric field strength as the boundary condtition. Then the Neumann boundary $\frac{\partial \phi}{\partial n} = E_{field}$ can be used.
 #. (ground): Dirichlet boundary $\phi = 0$.
 #. (insulation): Neumann boundary $\frac{\partial \phi}{\partial n} = 0$.

Weak Form of the Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^

To implement the :eq:`nernstplanck` and :eq:`poisson` in Hermes2D, the weak form must be derived. First of all let's denote:

* $K=z \mu F$
* $L=\frac{F}{\varepsilon}$

So equations :eq:`nernstplanck` and :eq:`poisson` can be written:

.. math::
	:label: nernstplancksimple
		
		\frac{\partial C}{\partial t}-D\Delta C-K\nabla\cdot \left(C\nabla\phi\right)=0,

.. math::
	:label: poissonsimple

		-\Delta\phi=L\left(C-C_{0}\right),

Then the boundary condition :eq:`nernstboundary` becomes

.. math::
	:label: nernstboundarysimple

		-D\frac{\partial C}{\partial n}-KC\frac{\partial\phi}{\partial n}=0.

Weak form of equation :eq:`nernstplancksimple` is:

.. math::
	:label: nernstweak1

		\int_{\Omega}\frac{\partial C}{\partial t}v d\mathbf{x}
		-\int_{\Omega}D\Delta Cv d\mathbf{x}-\int_{\Omega}K\nabla C\cdot
		\nabla\phi v d\mathbf{x} - \int_{\Omega}KC\Delta \phi v d\mathbf{x}=0,

where $v$ is a test function  $\Omega\subset\mathbf{R}^{3}$. When applying
Green's first identity to expand the terms that contain Laplacian
and adding the boundary condition :eq:`nernstboundarysimple`, the :eq:`nernstweak1`
becomes:

.. math::
	:label: nernstweak2

		\int_{\Omega}\frac{\partial C}{\partial t}v d\mathbf{x}+
		D\int_{\Omega}\nabla C\cdot\nabla v d\mathbf{x}-
		K\int_{\Omega}\nabla C \cdot \nabla \phi v d\mathbf{x}+
		K\int_{\Omega}\nabla\left(Cv\right)\cdot \nabla \phi d\mathbf{x}-
		D\int_{\Gamma}\frac{\partial C}{\partial n}v d\mathbf{S}-
		\int_{\Gamma}K\frac{\partial\phi}{\partial n}Cv d\mathbf{S}=0,

where the terms 5 and 6 equal $0$ due to the boundary condition. 
By expanding the nonlinear 4th term, the weak form becomes:

.. math::
	:label: nernstweak3

		\int_{\Omega}\frac{\partial C}{\partial t}v d\mathbf{x}+
		D\int_{\Omega}\nabla C \cdot \nabla v d\mathbf{x}-
		K\int_{\Omega}\nabla C \cdot \nabla \phi v d\mathbf{x}+
		K\int_{\Omega}\nabla \phi \cdot \nabla C v d\mathbf{x}+
		K\int_{\Omega} C \left(\nabla\phi\cdot\nabla v\right) d\mathbf{x}=0

As the terms 3 and 4 are equal and cancel out, the final weak form of equation
:eq:`nernstplancksimple` is

.. math::
	:label: nernstweak4

		\int_{\Omega}\frac{\partial C}{\partial t}v d\mathbf{x}+
		D\int_{\Omega}\nabla C \cdot \nabla v d\mathbf{x}+
		K\int_{\Omega} C \left(\nabla\phi\cdot\nabla v\right) d\mathbf{x}=0
		
The weak form of equation :eq:`poissonsimple` with test function $u$ is:

.. math::
	:label: poissonweak1

		-\int_{\Omega}\Delta\phi u d\mathbf{x}-\int_{\Omega}LCu d\mathbf{x}+
		\int_{\Omega}LC_{0}u d\mathbf{x}=0.

After expanding the Laplace' terms, the equation becomes:

.. math::
	:label: poissonweak2

		\int_{\Omega}\nabla\phi\cdot\nabla u d\mathbf{x}-\int_{\Omega}LCu d\mathbf{x}+
		\int_{\Omega}LC_{0}u d\mathbf{x}=0.

Notice, when electric field strength is used as a boundary condition, then the contribution of
the corresponding surface integral must be added:

.. math::
	:label: poissonweak3

		\int_{\Omega}\nabla\phi\cdot\nabla u d\mathbf{x}-\int_{\Omega}LCu d\mathbf{x}+
		\int_{\Omega}LC_{0}u d\mathbf{x}+\int_{\Gamma}\frac{\partial \phi}{\partial n}u d\mathbf{S}=0.

However, for the most cases we use only Poisson boundary conditions to set the voltage. Therefore the last
term of :eq:`poissonweak3` is omitted and :eq:`poissonweak2` is used instead in the following sections.

Jacobian matrix
^^^^^^^^^^^^^^^

Equation :eq:`nernstweak3` is time dependent, thus some time stepping 
method must be chosen. For simplicity we start with first order Euler implicit method

.. math::
	:label: euler

		\frac{\partial C}{\partial t} \approx \frac{C^{n+1} - C^n}{\tau}

where $\tau$ is the time step. We will use the following notation:

.. math::
	:label: cplus

		C^{n+1} = \sum_{k=1}^{N^C} y_k^{C} v_k^{C}, \ \ \ 
		  \phi^{n+1} = \sum_{k=1}^{N^{\phi}} y_k^{\phi} v_k^{\phi}.

In the new notation, time-discretized equation :eq:`nernstweak4` becomes:

.. math::
	:label: Fic

		F_i^C(Y) = \int_{\Omega} \frac{C^{n+1}}{\tau}v_i^C d\mathbf{x} - 
		\int_{\Omega} \frac{C^{n}}{\tau}v_i^C d\mathbf{x}
		+ D\int_{\Omega} \nabla C^{n+1} \cdot \nabla v_i^C d\mathbf{x}  
		+ K \int_{\Omega}C^{n+1} (\nabla \phi^{n+1} \cdot \nabla v_i^C) d\mathbf{x},

and equation :eq:`poissonweak2` becomes:

.. math::
	:label: Fiphi

		F_i^{\phi}(Y) = \int_{\Omega} \nabla \phi^{n+1} \cdot \nabla v_i^{\phi} d\mathbf{x} 
		- \int_{\Omega} LC^{n+1}v_i^{\phi} d\mathbf{x} + \int_{\Omega} LC_0 v_i^{\phi} d\mathbf{x}.

The Jacobian matrix $DF/DY$ has $2\times 2$ block structure, with blocks 
corresponding to

.. math:: 
	:label: jacobianelements

		\frac{\partial F_i^C}{\partial y_j^C}, \ \ \ \frac{\partial F_i^C}{\partial y_j^{\phi}}, \ \ \ 
		\frac{\partial F_i^{\phi}}{\partial y_j^C}, \ \ \ \frac{\partial F_i^{\phi}}{\partial y_j^{\phi}}.

Taking the derivatives of $F^C_i$ with respect to $y_j^C$ and $y_j^{\phi}$, we get

.. math::
	:label: bilin1

		\frac{\partial F_i^C}{\partial y_j^C} = 
		\int_{\Omega} \frac{1}{\tau} v_j^C v_i^C d\mathbf{x} + 
		D\int_{\Omega} \nabla v_j^C \cdot \nabla v_i^C d\mathbf{x}
		+ K\int_{\Omega} v_j^C (\nabla \phi^{n+1} \cdot \nabla v_i^C) d\mathbf{x},
	
.. math::
	:label: bilin2
		
		\frac{\partial F_i^C}{\partial y_j^{\phi}} =
		K \int_{\Omega} C^{n+1} (\nabla v_j^{\phi} \cdot \nabla v_i^C) d\mathbf{x}.

Taking the derivatives of $F^{\phi}_i$ with respect to $y_j^C$ and $y_j^{\phi}$, we get

.. math::
	:label: bilin3
		
		\frac{\partial F_i^{\phi}}{\partial y_j^C} =
		- \int_{\Omega} L v_j^C v_i^{\phi} d\mathbf{x},

.. math::
	:label: bilin4
		
		\frac{\partial F_i^{\phi}}{\partial y_j^{\phi}} =
		\int_{\Omega} \nabla v_j^{\phi} \cdot \nabla v_i^{\phi} d\mathbf{x}.

In Hermes, equations :eq:`Fic` and :eq:`Fiphi` are used to define the residuum $F$, and
equations :eq:`bilin1` - :eq:`bilin4` to define the Jacobian matrix $J$.
It must be noted that in addition to the implicit Euler iteration Crank-Nicolson iteration is implemented 
in the code (see the next section for the references of the source files).

Simulation
^^^^^^^^^^

To begin with simulations in Hermes2D, the equations :eq:`Fic` - :eq:`bilin4` were be implemented.
It was done by implementing the callback functions found in  `nernst-planck-timedep-adapt/forms.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/examples/nernst-planck/nernst-planck-timedep-adapt/forms.cpp>`_.

.. highlight:: c

The functions along with the boundary conditions::

	// Poisson takes Dirichlet and Neumann boundaries
	BCType phi_bc_types(int marker) {
		  return (marker == SIDE_MARKER || (marker == TOP_MARKER && VOLT_BOUNDARY == 2))
		      ? BC_NATURAL : BC_ESSENTIAL;
	}

	// Nernst-Planck takes Neumann boundaries
	BCType C_bc_types(int marker) {
		  return BC_NATURAL;
	}

	// Diricleht Boundary conditions for Poisson equation.
	scalar essential_bc_values(int ess_bdy_marker, double x, double y) {
		  return ess_bdy_marker == TOP_MARKER ? VOLTAGE : 0.0;
	}

are assembled as follows::
	
        // Add the bilinear and linear forms.
        if (TIME_DISCR == 1) {  // Implicit Euler.
          wf.add_vector_form(0, callback(Fc_euler), H2D_ANY,
		             Tuple<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton));
          wf.add_vector_form(1, callback(Fphi_euler), H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
          wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), HERMES_NONSYM, H2D_ANY, &phi_prev_newton);
          wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), HERMES_NONSYM, H2D_ANY, &C_prev_newton);
          wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), HERMES_NONSYM);
          wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), HERMES_NONSYM);
        } else {
          wf.add_vector_form(0, callback(Fc_cranic), H2D_ANY, 
		             Tuple<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton, &phi_prev_time));
          wf.add_vector_form(1, callback(Fphi_cranic), H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
          wf.add_matrix_form(0, 0, callback(J_cranic_DFcDYc), HERMES_NONSYM, H2D_ANY, Tuple<MeshFunction*>(&phi_prev_newton, &phi_prev_time));
          wf.add_matrix_form(0, 1, callback(J_cranic_DFcDYphi), HERMES_NONSYM, H2D_ANY, Tuple<MeshFunction*>(&C_prev_newton, &C_prev_time));
          wf.add_matrix_form(1, 0, callback(J_cranic_DFphiDYc), HERMES_NONSYM);
          wf.add_matrix_form(1, 1, callback(J_cranic_DFphiDYphi), HERMES_NONSYM);
        }

where the variables ``C_prev_time``, ``C_prev_newton``, 
``phi_prev_time``, and ``phi_prev_newton`` are solutions of concentration
$C$ and voltage $\phi$. The suffixes *newton* and *time* are current iteration and previous
time step, respectively.

When it comes to meshing, it should be considered that the gradient of $C$ near the boundaries will
be higher than gradients of $\phi$. This allows us to create different meshes for those variables. In
`main.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/examples/nernst-planck/nernst-planck-timedep-adapt/main.cpp>`_.
the following code in the *main()* function enables multimeshing


.. code-block:: c
	
    // Spaces for concentration and the voltage.
    H1Space C(&Cmesh, C_bc_types, C_essential_bc_values, P_INIT);
    H1Space phi(MULTIMESH ? &phimesh : &Cmesh, phi_bc_types, phi_essential_bc_values, P_INIT);

When ``MULTIMESH`` is defined in `main.cpp <http://git.hpfem.org/hermes.git/blob/HEAD:/hermes2d/examples/nernst-planck/nernst-planck-timedep-adapt/main.cpp>`_.
then different H1Spaces for ``phi`` and ``C`` are created. It must be noted that when adaptivity
is not used, the multimeshing in this example does not have any advantage, however, when
adaptivity is turned on, then mesh for H1Space ``C`` is refined much more than for ``phi``.

Non adaptive solution
^^^^^^^^^^^^^^^^^^^^^

The following figure shows the calculated concentration $C$ inside the IPMC.

.. image:: example-np/nonadapt_conc.png
	:align: center
	:alt: Calculated concentration

As it can be seen, the concentration is rather uniform in the middle of domain. In fact, most of the
concentration gradient is near the electrodes, within 5...10% of the total thickness. That is why the refinement

The voltage inside the IPMC forms as follows:

.. image:: example-np/nonadapt_phi.png
	:align: center
	:alt: Calculated voltage inside the IPMC

Here we see that the voltage gradient is smaller and more uniform near the boundaries than it is for $C$.
That is where **the adaptive multimeshing** can become useful.

Adaptive solution
^^^^^^^^^^^^^^^^^

To be added soon.

