Essential and Natural Boundary Conditions
-----------------------------------------

Hermes distinguishes between **essential** and **natural** boundary conditions. 
The former type eliminates degrees of freedom from the domain's boundary
(the solution or its portion is prescribed) while the latter does not. 
Examples of essential boundary conditions are Dirichlet conditions for 
$H^1$ problems and perfect conductor conditions in the space $H(curl)$.
Examples of natural boundary conditions are Neumann or Newton (Robin) 
conditions for $H^1$ problems and impedance conditions in the space 
$H(curl)$. Only essential conditions are treated explicitly in Hermes, 
while the natural ones the user builds into the surface integrals 
of his/her weak formulation. Let us start with showing default ways 
to define essential boundary conditions.

Default constant essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We already met a default essential boundary condition in the previous example 
`03-poisson <http://hpfem.org/hermes/doc/src/hermes2d/P01-linear/03-poisson.html>`_::

    // Initialize essential boundary conditions.
    DefaultEssentialBCConst bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), FIXED_BDY_TEMP);

This one assigned a constant value FIXED_BDY_TEMP to all boundary edges with the markers 
"Bottom", "Inner", "Outer" or "Left". 

After creating one or more essential boundary conditions, they are passed into the container 
class EssentialBCs::

    EssentialBCs bcs(&bc_essential);

The purpose of this container is to collect all essential boundary conditions to be passed into a Space, 
as we shall see shortly. The constructor of the EssentialBCs class can accept a Hermes::vector of
essential boundary conditions. 

Default nonconstant essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another default essential boundary condition reads boundary values from a given 
function. This is useful, for example, in benchmarks with known exact solution. A typical
usage is as follows (taken from `benchmark NIST-01 <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks-nist/nist-01.html>`_)::

    // Set exact solution.
    CustomExactSolution exact(&mesh, EXACT_SOL_P);

    // Initialize boundary conditions
    DefaultEssentialBCNonConst bc_essential("Bdy", &exact);
    EssentialBCs bcs(&bc_essential);

This will be discussed in more detail in the benchmarks section.

Custom essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Custom essential conditions can be created by subclassing the abstract class
EssentialBoundaryCondition::

    class HERMES_API EssentialBoundaryCondition
    {
    public:
      /// Default constructor.
      EssentialBoundaryCondition(Hermes::vector<std::string> markers);
      EssentialBoundaryCondition(std::string marker);

      /// Virtual destructor.
      virtual ~EssentialBoundaryCondition();

      /// Types of description of boundary values, either a function (callback), or a constant.
      enum EssentialBCValueType {
	BC_FUNCTION,
	BC_CONST
      };

      /// Pure virtual function reporting the type of the essential boundary condition.
      virtual EssentialBCValueType get_value_type() const = 0;

      /// Represents a function prescribed on the boundary. Gets the boundary point coordinate as well as the 
      /// normal and tangential vectors.
      virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const = 0;

      /// Special case of a constant function.
      scalar value_const;

      /// Sets the current time for time-dependent boundary conditions.
      void set_current_time(double time);
      double get_current_time() const;

    protected:
      /// Current time.
      double current_time;

      // Markers.
      Hermes::vector<std::string> markers;

      // Friend class.
      friend class EssentialBCs;
      friend class Space;
    };

This class can represent arbitrary essential boundary conditions that depend 
on space and time. Every descendant of this class must redefine the purely 
virtual functions get_value_type() and value(). This will be explained in
more detail in the following example.  

