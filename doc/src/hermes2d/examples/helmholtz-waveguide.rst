Helmholtz equation - waveguides (Electromagnetics)
------------------------------

Mathematical description of waveguides
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The most general mathematical description of the phenomena in waveguides is given by Maxwell equations  

.. math::
	:label: 1. Maxwell equation
		
		\nabla \times {\pmb{H}} &= {\pmb{J}} +
		\frac{\partial {\pmb{D}}}{\partial t}, 

.. math::
	:label: 2. Maxwell equation	
		
		\nabla \times {\pmb{E}} &= 
		- \frac{\partial {\pmb{B}}}{\partial t},
	
.. math::
	:label: 3. Maxwell equation		
		
		\nabla \cdot {\pmb{D}} &= \rho, 
		
.. math::
	:label: 4. Maxwell equation		
		
		\nabla \cdot {\pmb{B}} &= 0, 	

with respect

.. math::
	:label: Material properties
	
		{\pmb{B}} = \mu {\pmb{H}}, \ 
		{\pmb{J}} = \sigma {\pmb{E}}, \
		{\pmb{D}} = \varepsilon {\pmb{E}},
		
where  :math:`\varepsilon` means permitivity, :math:`\mu` permeability and :math:`\sigma` stands for electric conductivity. For the waveguides analysis purpose material properties are often considered as constants and isotropic. After substituting material properties :eq:`Material properties` into equations :eq:`1. Maxwell equation` and :eq:`2. Maxwell equation`,  we get

.. math::
	:label: 1.a Maxwell equation	

		 \nabla \times \frac{1}{\mu} {\pmb{B}} &= \sigma {\pmb{E}} +
		\varepsilon \frac{\partial {\pmb{E}}}{\partial t}, 

.. math::
	:label: 2.a Maxwell equation	

		 \nabla \times {\pmb{E}} &= 
		- \frac{\partial {\pmb{B}}}{\partial t}. 

If the vector operator :math:`\mathrm{curl}` is applied on the equation :eq:`2.a Maxwell equation`, it is possible to substitute :math:`\nabla \times \pmb{E}` from the equation :eq:`1.a Maxwell equation` and get the wave equation for the electric field in the form

.. math::
	:label: Wave equation
	
		\nabla \times \nabla \times \pmb{E} =
		- \mu \sigma \frac{\partial {\pmb{E}}}{\partial t} 
		- \mu \varepsilon \frac{\partial^2 {\pmb{E}}}{\partial t^2}. 

In a medium with the zero charge density :math:`\rho` it is usefull to apply the vector identity (because :math:`\nabla \cdot \pmb{E} = 0`)

.. math::
	:label: 1. vector identity
	
		\nabla \times \nabla \times \pmb{E} = \nabla \nabla \cdot \pmb{E} - \Delta \pmb{E}


The wave equation :eq:`Wave equation` can be siplified to the most common form

.. math::
	:label: a. Wave equation
	
		\Delta \pmb{E} - \mu \sigma \frac{\partial {\pmb{E}}}{\partial t} - \mu \varepsilon \frac{\partial^2 {\pmb{E}}}{\partial t^2} = \mathbf{0}.
	
For many technical problems it is sufficient to know the solution in the frequency domain. After applying Fourier transform, the equation :eq:`a. Wave equation` becomes 

.. math::
	:label: Hemlholtz equation

	\Delta \overline{\pmb{E}} - \mathrm{j} \mu \sigma \omega \overline{\pmb{E}} + \omega^2 \mu \varepsilon \overline{{\pmb{E}}} = \mathbf{0},

	
which is in fact the Hemlholtz equation.

Parallel plate waveguide - harmonic analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The parallel plate waveguide is the simplest type of guide that supports TM (transversal magnetic) and TE (transversal electric) modes. This kind of guide allows also TEM (transversal elektric and magnetic) mode.

Geometry
^^^^^^^^

.. image:: hemlholtz-waveguide/waveguide.png
   :scale: 50 %   
   :align: center 	
   :alt: Paralel plate waveguide geometry
	
Mathematical model - TE modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose that the electromagnetic wave is propagating in the direction :math:`z`, then the component of the vector :math:`\pmb{E}` in the direction of the propagation is equal to zero

.. math::
	:label: a. TE mode

	\overline{E_z} = 0,

thus it is possible to solve the electric field in the parallel plate waveguide as a two-dimensional problem. 

On the whole area the wave equation (Helmholtz equation) is fulfilled 

.. math::
    :label: a. Helmholtz equation

       \Delta \overline{\pmb{E}} - \mathrm{j} \mu \sigma \omega \overline{\pmb{E}} + \omega^2 \mu \varepsilon \overline{{\pmb{E}}} = \mathbf{0}.

The conducting plates (boundary :math:`\Gamma_1, \Gamma_2`) are usually supposed as *perfectly conductive*, that can be modelled using the boundary condition

.. math::
	:label: Perfect conductor

	\pmb{n} \times \overline{\pmb{E}} = 0.

for the geometry on the picture the whole expression :eq:`Perfect conductor` is reduced to the zero Dirichlet boundary condition

.. math::
		:label: Reduced Perfect conductor

		\overline{E_x} = 0.


For the bounderies :math:`\Gamma_3, \Gamma_4` there are commonly used several types of boundery conditions:

Electric field (Dirichlet boundary condition)
"""""""""""""""""""""""""""""""""""""""""""""

	.. math::
		:label: Electric field

			\overline{\pmb{E}}(\Gamma) = \overline{E_0} = \mathrm{const}.

	Note that for TE modes (and for the geometry on the picture) the natural boundery condition is described by the expression

	.. math::
		:label: TE Electric field

		\overline{E}_x(y) = \overline{E_0} \cos\left(\frac{y \cdot n \pi}{h} \right),

	where :math:`n` stands for mode.

Impedance matching (Newton boundary condition)
""""""""""""""""""""""""""""""""""""""""""""""

	For harmonic TE mode waves the following relation is valid

	.. math::
		:label: Impedance definition

		\overline{\pmb{E}} = Z_0 (\overline{H_y} \pmb{i} - \overline{H_x} \pmb{j}) = Z_0 \cdot \pmb{n} \times \overline{\pmb{H}},

	where :math:`Z_0` is *the wave impedance*. At the same time the second Maxwell equation

	.. math::
		:label: Harmonic Maxwell equation

		\nabla\times \overline{{\pmb{E}}} = -j \omega \mu \overline{\pmb{H}}
	
	must be satisfied. From quations :eq:`Impedance definition` and :eq:`Harmonic Maxwell equation` it is possible to derive impedance matching boundary condition in the form

	.. math::
		:label: Impedance matching

		\pmb{n} \times \nabla \times \overline{\pmb{E}} =  \frac{j \omega \mu }{Z_0} \overline{\pmb{E}} =  j \beta \overline{\pmb{E}}.

	For given geometry the equation :eq:`Impedance matching` can be reduced to the Newton boundary condition in a form

	..  math::
		:label: Newton boundary condition

		\frac{\partial \overline{E_x}}{\partial y} = j \beta \overline{E_x}.

Program
^^^^^^^
Material parameters
"""""""""""""""""""
::

	const double epsr = 1.0;                    // Relative permittivity
	const double eps0 = 8.85418782e-12;         // Permittivity of vacuum F/m
	const double mur = 1.0;                     // Relative permeablity
	const double mu0 = 4*M_PI*1e-7;             // Permeability of vacuum H/m
	const double frequency = 3e9;               // Frequency MHz
	const double omega = 2*M_PI * frequency;    // Angular velocity
	const double sigma = 0;                     // Conductivity Ohm/m

Boundary conditions
"""""""""""""""""""
There are three types of boundary conditions:

	* Zero Dirichlet boundary condition on boundaries :math:`\Gamma_1, \Gamma_2`:
	::

		bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_PERFECT, BDY_LEFT));
		BCValues bc_values_r;
		bc_values_r.add_const(BDY_PERFECT, 0.0);
		BCValues bc_values_i;
		bc_values_i.add_const(BDY_PERFECT, 0.0);

	Note: There are two ways how to approach to complex problems. The first one is using the type complex from the STL library. The second approach is demonstrated - two matrices and two right side vectors (one for real parts and second for imagnary parts) are build.

	* Nonzero Dirichlet boundary conditions:
		
	The boundary condition on the boundary :math:`\Gamma_3` is specified according to the expression :eq:`TE Electric field`.::

		scalar essential_bc_values(double x, double y)
		{
		  return cos(y*M_PI/0.1)*100;
		}

		int main()
		{
		  ...
		  bc_values_r.add_function(BDY_LEFT, essential_bc_values);
		  ...		  
		}

	* Newton boundary condition
	::

		bc_types.add_bc_newton(Hermes::vector<int>(BDY_IMPEDANCE));


Weak forms
""""""""""
	* registration (in function ``main()``) ::

		WeakForm wf(2);
		wf.add_matrix_form(0, 0, callback(magnetic_matrix_form_real_real));
		wf.add_matrix_form(0, 1, callback(magnetic_matrix_form_real_imag));
		wf.add_matrix_form(1, 1, callback(magnetic_matrix_form_imag_imag));
		wf.add_matrix_form(1, 0, callback(magnetic_matrix_form_imag_real));
		wf.add_matrix_form_surf(0, 1, callback(magnetic_vector_form_surface_imag_real), BDY_IMPEDANCE);
		wf.add_matrix_form_surf(1, 0, callback(magnetic_vector_form_surface_real_imag), BDY_IMPEDANCE);
	
	The function ``magnetic_matrix_form_real_real`` describes behaviour of the of the component of electric field :math:`\overline{E_x}` and ``magnetic_matrix_form_imag_imag`` describes behaviour of the imaginary part of the component of electric field :math:`\overline{E_x}` in this code. Functions ``magnetic_matrix_form_imag_real`` and ``magnetic_matrix_form_real_imag`` represent the conection between real and imaginary part. The functions ``magnetic_vector_form_surface_imag_real`` and ``magnetic_vector_form_surface_imag_real`` express the Newton boundary condition and also the conection between the real and the imaginary part of the component of the electric field :math:`\overline{E_x}`. 

Results
^^^^^^^
.. image:: helmholtz-waveguide/real_part.png
   :scale: 50 %   
   :align: center 	
   :alt: Paralel plate waveguide geometry

.. image:: helmholtz-waveguide/imaginary_part.png
   :scale: 50 %   
   :align: center 	
   :alt: Paralel plate waveguide geometry
	
