Tutorial Part IV (Adaptive Solution of Linear Problems)
=======================================================

So far we have not paid any attention to the accuracy of results. In general, 
a computation on a fixed mesh is not likely to be very accurate. There is a need 
for *adaptive mesh refinement (AMR)* that improves the quality of the approximation 
by refining mesh elements or increases the polynomial degree where the approximation 
error is large.

.. toctree::
    :maxdepth: 2

    linear-adapt/adapt   
    linear-adapt/conv   
    linear-adapt/micromotor 
    linear-adapt/multimesh
    linear-adapt/multimesh-example
    linear-adapt/general-adapt
    linear-adapt/complex-adapt
    linear-adapt/maxwell-adapt
    linear-adapt/adapt_exact

 








