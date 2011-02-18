Newton's Method
---------------

The Newton's method is more powerful but also more 
demanding than the Picard's method, and therefore 
we will begin at a slower pace. 

We'll stay with the model problem from the previous section

.. math::
    :label: newton0

    -\nabla \cdot (\lambda(u)\nabla u) - f(\bfx) = 0, \ \ \ u = 0 \ \mbox{on}\ \partial \Omega.

Note that when using the Newton's method, it is customary to have 
everything on the left-hand side. The corresponding discrete problem has the form 

.. math::

    \int_{\Omega} \lambda(u)\nabla u(\bfx) \cdot \nabla v_i(\bfx)\, \mbox{d}\bfx 
    - \int_{\Omega} f(\bfx)v_i(\bfx) \, \mbox{d}\bfx = 0\ \ \ \mbox{for all} \ i = 1, 2, \ldots, N, 

where $v_i$ are the test functions and

.. math::

    u(\bfY) = \sum_{j=1}^N y_j v_j.

Here $\bfY = (y_1, y_2, \ldots, y_N)^T$ is the vector of unknown coefficients.
The nonlinear discrete problem can be written in the compact form

.. math::

    \bfF(\bfY) = {\bf 0},
 
where $\bfF = (F_1, F_2, \ldots, F_N)^T$ is the so-called *residuum* or *residual vector*.

Residual vector and Jacobian matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The residual vector $\bfF(\bfY)$ has $N$ components of the form

.. math::

    F_i(\bfY) =  \int_{\Omega} \lambda(u)\nabla u \cdot \nabla v_i 
    - f v_i \, \mbox{d}\bfx.

Here $v_i$ is the $i$-th test function, $i = 1, 2, \ldots, N$.
The $N\times N$ Jacobian matrix $\bfJ(\bfY) = D\bfF/D\bfY$ has the components 

.. math::

    J_{ij}(\bfY) =  \frac{\partial F_i}{\partial y_j} = 
    \int_{\Omega} \left[ \frac{\partial \lambda}{\partial u} \frac{\partial u}{\partial y_j} 
    \nabla u + \lambda(u)\frac{\partial \nabla u}{\partial y_j} \right] \cdot \nabla v_i \, \mbox{d}\bfx.

Using elementary relations shown below, we obtain

.. math::

    J_{ij}(\bfY) =
    \int_{\Omega} \left[ \frac{\partial \lambda}{\partial u}(u) v_j 
    \nabla u + \lambda(u)\nabla v_j \right] \cdot \nabla v_i \, \mbox{d}\bfx.

It is worth noticing that $\bfJ(\bfY)$ has the same sparsity structure as the 
standard stiffness matrix that we know from linear problems. In fact, when the 
problem is linear then the Jacobian matrix and the stiffness matrix are the same 
thing (see last paragraph in this section). 

How to differentiate weak forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chain rule is used all the time. If your equation contains parameters that depend on 
the solution, you will need their derivatives with respect to the solution (such as we needed 
the derivative of $\lambda$ above). In addition, the following elementary rules are useful 
for the differentiation of the weak forms: 

.. math::

    \frac{\partial u}{\partial y_k} = \frac{\partial}{\partial y_k}\sum_{j=1}^N y_j v_j = v_k

and 

.. math::

    \frac{\partial \nabla u}{\partial y_k} = \frac{\partial}{\partial y_k}\sum_{j=1}^N y_j \nabla v_j = \nabla v_k.

Practical formula
~~~~~~~~~~~~~~~~~

Let's assume that the Jacobian matrix has been assembled. 
The Newton's method is written formally as 

.. math::

    \bfY_{\!\!n+1} = \bfY_{\!\!n} - \bfJ^{-1}(\bfY_{\!\!n}) \bfF(\bfY_{\!\!n}),

but a more practical formula to work with is 

.. math::

    \bfJ(\bfY_{\!\!n})\delta \bfY_{\!\!n+1} =  - \bfF(\bfY_{\!\!n}).

This is a system of linear algebraic equations that needs to be solved in every Newton's 
iteration. The Newton's method will stop when $\bfF(\bfY_{\!\!n+1})$ is sufficiently close 
to the zero vector.

Stopping criteria
~~~~~~~~~~~~~~~~~

There are two basic stopping criteria that should be combined 
for safety:

* Checking whether the (Euclidean) norm of the residual vector $\bfF(\bfY_{n+1})$ is sufficiently close to zero.
* Checking whether the norm of $\bfY_{n+1} - \bfY_{n}$ is sufficiently close to zero.

If just one of these two criteria is used, the Newton's method may finish prematurely.

Linear problems
~~~~~~~~~~~~~~~

In the linear case we have 

.. math::

    \bfF(\bfY) = \bfJ(\bfY)\bfY - \bfb,

where $\bfS = \bfJ(\bfY)$ is a constant stiffness matrix and $\bfb$ a load vector. 
The Newton's method is now

.. math::

    \bfS\bfY_{\!\!n+1} = \bfJ(\bfY_{\!\!n})\bfY_{\!\!n} 
    - \bfJ(\bfY_{\!\!n})\bfY_{\!\!n} + \bfb = \bfb.

There exists a widely adopted mistake saying that 
the Newton's method, when applied to a linear problem, 
will converge in one iteration. This is only true if 
one uses the first (residual norm based) stopping 
criterion above. If the second criterion is used, 
which is based on the distance of two consecutive 
solution vectors, then the Newton's method will do 
two steps before stopping. In practice, using just 
the residual criterion is dangerous.

This explains that it makes sense to 
use the knowledge that the problem is linear, and 
stop the Newton's iteration after the first step 
manually.


