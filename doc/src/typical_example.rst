Hermes typical example structure
--------------------------------

- a nice thing about Hermes is that it follows the math tightly, so the steps taken in solving the example with Hermes correspond to those taken in theory.

A beginning of each example can look like this::
    
    // Include the main Hermes2D header.
    include "hermes2d.h"
    
    // Two basic namespaces.
    using namespace Hermes;
    using namespace Hermes::Hermes2D;
    
    // For adaptivity.
    using namespace Hermes::Hermes2D::RefinementSelectors;
    
    // For visualization.
    using namespace Hermes::Hermes2D::Views;

.. toctree::
    :maxdepth: 2

    typical_example/mesh
    typical_example/space
    typical_example/weak_formulation
    typical_example/calculation
    