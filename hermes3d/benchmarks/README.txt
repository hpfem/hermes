Benchmarks are slow!
--------------------
Most benchmarks come with a non-constant non-polynomial 
right hand side (enforced by an exact solution), and thus 
Hermes uses numerical quadrature of high orders to integrate 
it with sufficient accuracy. For testing purposes, the order 
of the quadrature can be always lowered manually in the 
corresponding weak forms. 

