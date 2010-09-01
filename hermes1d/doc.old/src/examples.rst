========
Examples
========

The best way to understand the above machinery is to solve examples which we
will do in this section.

Classical Harmonic Oscillator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


One of the important equations from the classical mechanics is the harmonic
oscillator equation:

.. math::

    u''(x)+u(x)=0

and for this example we choose a simple boundary conditions $u(0)=0$ and
$u'(0)=1$ so the solution is $u(x)=\sin
x$. First let's rewrite the equation into the form :eq:`one`:

.. math::

     {\d u_1\over\d x}=g_1(u_1, u_2, x)=u_2


.. math::

     {\d u_2\over\d x}=g_2(u_1, u_2, x)=-u_1

where $u_1=u$ is the function we seek and $u_2=u'$ is its derivative.
Then let's write it in the form :eq:`two`:

.. math::

     f_1(u_1, u_2, u'_1, u'_2, x) = g_1(u_1, u_2, x)-u_1'=-u_1'+u_2 = 0,


.. math::

     f_2(u_1, u_2, u'_1, u'_2, x) = g_2(u_1, u_2, x)-u_2'=-u_2'-u_1 = 0,

and :eq:`qqq`:

.. math::

     \bfF(\bfU, \bfU', x) = 0

where ${\bfU} = (u_1, u_2)$ and ${\bfF} = (f_1, f_2)=(-u_1'+u_2, -u_2'-u_1)$.
The weak formulation is:

.. math::

     R_1(\bfY) = \int_0^T f_1(u_1, u_2, u'_1, u'_2, x)v_1 \, \d t = \int_0^T (-u_1'+u_2)v_1 \, \d t =0,


.. math::

     R_2(\bfY) = \int_0^T f_2(u_1, u_2, u'_1, u'_2, x)v_1 \, \d t = \int_0^T (-u_2'-u_1)v_2 \, \d t =0,

To evaluate the Jacobi matrix $D\bfR/D\bfY$ for the Newton's iteration, we need
the following Jacobians:

.. math::

     \left({{\rm D}\bfF\over{\rm D}\bfU}\right)_{mn}= \frac{\partial f_m}{\partial u_n}(u_1, u_2, u'_1, u'_2, x) = \left( \begin{array}{c} -u_1'+u2 \\ -u_2'-u1 \\ \end{array} \right) \left( \begin{array}{cc} \overleftarrow{\partial_{u_1}} & \overleftarrow{\partial_{u_2}} \\ \end{array} \right) = \left( \begin{array}{cc} 0 & 1 \\ -1 & 0 \\ \end{array} \right)


.. math::

     \left({{\rm D}\bfF\over{\rm D}\bfU'}\right)_{mn}= \frac{\partial f_m}{\partial u'_n}(u_1, u_2, u'_1, u'_2, x) = \left( \begin{array}{c} -u_1'+u2 \\ -u_2'-u1 \\ \end{array} \right) \left( \begin{array}{cc} \overleftarrow{\partial_{u_1'}} & \overleftarrow{\partial_{u_2'}} \\ \end{array} \right) = \left( \begin{array}{cc} -1 & 0 \\ 0 & -1 \\ \end{array} \right)

where $\overleftarrow{\partial_{u_1}}$ is a partial derivative with respect to
$u_1$ but acting to the left.

To solve this problem with Hermes, all we have to do is to specify the
following information:

.. math::

    {\bfF} = \left( \begin{array}{c} -u_1'+u2 \\ -u_2'-u1 \\ \end{array} \right)


.. math::

     {{\rm D}\bfF\over{\rm D}\bfU}= \left( \begin{array}{cc} 0 & 1 \\ -1 & 0 \\ \end{array} \right)


.. math::

     {{\rm D}\bfF\over{\rm D}\bfU'}= \left( \begin{array}{cc} -1 & 0 \\ 0 & -1 \\ \end{array} \right)


.. math::

    u_1(0)=0


.. math::

    u_2(0)=1

and Hermes will solve for $\bfF=0$. This is implemented in
``examples/sin.py``.

Quantum Harmonic Oscillator
~~~~~~~~~~~~~~~~~~~~~~~~~~~


The corresponding quantum mechanics problem to the previous one is the quantum
harmonic oscillator for one particle in 1D:

.. math::

     i\hbar{\partial\over\partial t}\psi(x, t)= -{\hbar^2\over2m}{\partial^2\over\partial x^2}\psi(x,t)+V(x)\psi(x,t)


.. math::

     V(x)={1\over2}m\omega^2x^2

This is a partial differential equation for the time evolution of the wave
function $\psi(x, t)$, but one method to solve it is the
eigenvalues expansion:

.. math::

    \psi(x,t) = \sum_E c_E\psi_E(x)e^{-{i\over\hbar}Et}

where the sum goes over the whole spectrum (for continuous spectrum the sum
turns into an integral), the $c_E$ coefficients are determined from the initial condition
and $\psi_E(x)$ satisfies the one dimensional one particle time independent
Schroedinger equation:

.. math::

     -{\hbar^2\over2m}{\d^2\over\d x^2}\psi_E(x)+V(x)\psi_E(x)=E\psi_E(x)

and this is just an ODE and thus can be solved with Hermes1D. There can be many
types of boundary conditions for this equation, depending on the physical
problem, but in our case we simply have $\lim_{x\to\pm\infty}\psi_E(x)=0$ and
the normalization condition $\int_{-\infty}^\infty|\psi_E(x)|^2\d x=1$.

We can set $m=\hbar=1$ and from now on we'll just write $\psi(x)$ instead of
$\psi_E(x)$:

.. math::

     -{1\over2}{\d^2\over\d x^2}\psi(x)+V(x)\psi(x)=E\psi(x)

and we will solve it on the interval $(a, b)$ with the boundary condition
$\psi(a)=\psi(b)=0$. The weak formulation is

.. math::

     \int_a^b{1\over2}{\d\psi(x)\over\d x}{\d v(x)\over\d x}+V(x)\psi(x)v(x)\,\d x -\left[{\d\psi(x)\over\d x}v(x)\right]^a_b =E\int_a^b\psi(x)v(x)\,\d x

but due to the boundary condition $v(a)=v(b)=0$ so
$\left[\psi'(x)v(x)\right]^a_b=0$ and we get

.. math::

     \int_a^b{1\over2}{\d\psi(x)\over\d x}{\d v(x)\over\d x}+V(x)\psi(x)v(x)\,\d x =E\int_a^b\psi(x)v(x)\,\d x

And the finite element formulation is then $\psi(x)=\sum_j y_j\phi_j(x)$ and
$v=\phi_i(x)$:

.. math::

     \left(\int_a^b{1\over2}\phi_i'(x)\phi_j'(x)+V(x)\phi_i(x)\phi_j(x)\,\d x\right) y_j =E\int_a^b\phi_i(x)\phi_j(x)\,\d x\ y_j

which is a generalized eigenvalue problem:

.. math::

     A_{ij}y_j=EB_{ij}y_j

with

.. math::

     A_{ij}=\int_a^b{1\over2}\phi_i'(x)\phi_j'(x)+V(x)\phi_i(x)\phi_j(x)\,\d x


.. math::

     B_{ij}=\int_a^b\phi_i(x)\phi_j(x)\,\d x



Radial Schroedinger Equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Another important example is the three dimensional one particle time
independent Schroedinger equation for a spherically symmetric potential:

.. math::

     -{1\over2}\nabla^2\psi({\bf x})+V(r)\psi({\bf x})=E\psi({\bf x})

The way to solve it is to separate the equation into radial and angular parts
by writing the Laplace operator in spherical coordinates as:

.. math::

     \nabla^2f =  {\partial^2 f\over\partial\rho^2} +{2\over \rho}{\partial^2 f\over\partial\rho^2} -{L^2\over \rho^2}


.. math::

     L^2= -{\partial^2 f\over\partial\theta^2} -{1\over\sin^2\theta}{\partial^2 f\over\partial\phi^2} -{1\over\tan\theta}{\partial f\over\partial\theta}

Substituting $\psi=R(\rho)Y(\theta,\phi)$ into the Schroedinger equation
yields:

.. math::

    -{1\over2}\nabla^2(RY)+VRY=ERY


.. math::

    -{1\over2}R''Y-{1\over\rho}R'Y+{L^2RY\over2\rho^2}+VRY=ERY

Using the fact that $L^2Y=l(l+1)Y$ we can cancel $Y$ and we get the radial
Schroedinger equation:

.. math::

    -{1\over2}R''-{1\over\rho}R'+{l(l+1)R\over2\rho^2}+VR=ER

The solution is then:

.. math::

    \psi({\bf x})=\sum_{nlm}c_{nlm}R_{nl}(r)Y_{lm}\left({\bf x}\over r\right)

where $R_{nl}(r)$ satisfies the radial Schroedinger equation (from now on we
just write $R(r)$):

.. math::

    -{1\over2}R''(r)-{1\over r}R'(r)+\left(V+{l(l+1)\over2r^2}\right)R(r)=ER(r)

Again there are many types of boundary conditions, but the most common case is
$\lim_{r\to\infty}R(r)=0$ and $R(0)=1$ or $R(0)=0$. One solves this equation on
the interval $(0, a)$ for large enough $a$.

The procedure is similar to the previous example, only we need to remember that
we always have to use covariant integration (in the previous example the
covariant integration was the same as the coordinate integration),
in this case $r^2\sin\theta \d
r\d\theta\d\phi$, so the weak formulation is:

.. math::

    \int \left(-{1\over2}R''(r)-{1\over r}R'(r)+\left(V+{l(l+1)\over2r^2}\right)R(r)\right)v(r)r^2\sin\theta \d r\d\theta\d\phi=


.. math::

     =\int ER(r) v(r)r^2\sin\theta \d r\d\theta\d\phi

Integrating over the angles gives $4\pi$ which we cancel out at both sides and
we get:

.. math::

    \int_0^a \left(-{1\over2}R''(r)-{1\over r}R'(r)+\left(V+{l(l+1)\over2r^2}\right)R(r)\right)v(r)r^2 \d r=


.. math::

     =E\int_0^a R(r) v(r)r^2 \d r

We apply per partes to the first two terms on the left hand side:

.. math::

    \int_0^a \left(-{1\over2}R''(r)-{1\over r}R'(r)\right)v(r)r^2 \d r =\int_0^a -{1\over2r^2}\left(r^2 R'(r)\right)'v(r)r^2 \d r=


.. math::

     =\int_0^a -{1\over2}\left(r^2 R'(r)\right)'v(r) \d r =\int_0^a {1\over2}r^2 R'(r)v'(r) \d r-{1\over2} [r^2R'(r)v(r)]_0^a=


.. math::

     =\int_0^a {1\over2} R'(r)v'(r) r^2\d r -{1\over2} a^2R'(a)v(a)

We used the fact that $\lim_{r\to0} r^2 R'(r) = 0$. If we also prescribe the
boundary condition $R'(a)=0$, then the boundary term vanishes completely. The
weak formulation is then:

.. math::

    \int_0^a {1\over2}R'(r)v'(r)r^2+ \left(V+{l(l+1)\over2r^2}\right)R(r)v(r)r^2\,\d r = E\int_0^aR(r)v(r)r^2\,\d r

or

.. math::

    \int_0^a {1\over2}R'(r)v'(r)r^2+ V(r)R(r)v(r)r^2+{l(l+1)\over2} R(r)v(r)\,\d r = E\int_0^aR(r)v(r)r^2\,\d r


Another (equivalent) approach is to write a weak formulation for
the 3D problem in cartesian coordinates:

.. math::

     \int_\Omega{1\over2}\nabla\psi({\bf x})\nabla v({\bf x})+V(r)\psi({\bf x})v({\bf x})\,\d^3 x =E\int_\Omega\psi({\bf x})v({\bf x})\,\d^3 x

and only then transform to spherical coordinates:

.. math::

     \int_0^{2\pi}\d\varphi\int_0^\pi\d\theta\int_0^a\d r \left({1\over2}\nabla\psi({\bf x})\nabla v({\bf x})+V(r)\psi({\bf x})v({\bf x})\right)r^2\sin\theta=


.. math::

     = E\int_0^{2\pi}\d\varphi\int_0^\pi\d\theta\int_0^a\d r\, \psi({\bf x})v({\bf x})r^2\sin\theta

The 3d eigenvectors $\psi({\bf x})$ however are not spherically symmetric.
Nevertheless we can still proceed by choosing our basis as

.. math::

    v_{ilm}({\bf x})=\phi_{il}(r)Y_{lm}(\theta, \varphi)

and seek our solution as

.. math::

    \psi({\bf x})=\sum_{jlm}y_{jlm}\phi_{jl}(r)Y_{lm}(\theta, \varphi)

Using the properties of spherical harmonics and the gradient:

.. math::

    \int Y_{lm} Y_{l'm'} \sin\theta\,\d\theta\,\d\varphi= \delta_{ll'}\delta_{mm'}


.. math::

    \int r^2\nabla Y_{lm} \nabla Y_{l'm'} \sin\theta\,\d\theta\,\d\varphi= l(l+1)\delta_{ll'}\delta_{mm'}


.. math::

    \nabla f = {\partial f\over \partial r}\boldsymbol{\hat r} + {1\over r} {\partial f\over\partial\theta}\boldsymbol{\hat\theta}+{1\over r\sin\theta} {\partial f\over\partial\phi}\boldsymbol{\hat\phi}

the weak formulation becomes:

.. math::

     \left(\int_0^a {1\over2}r^2\phi_{il}'(r)\phi_{jl}'(r)+ {1\over2}X+ {l(l+1)\over2}\phi_{il}(r)\phi_{jl}(r)+ r^2V(r)\phi_{il}(r)\phi_{jl}(r)\,\d r\right)y_{jlm}=


.. math::

     = E\int_0^ar^2 \phi_{il}(r)\phi_{jl}(r)\,\d r\ y_{jlm}

where both $l$ and $m$ indices are given by the indices of the particular base
function $v_{ilm}$. The $X$ term is (schematically):

.. math::

    X=\int r^2\sin\theta(r)Y_{lm}(\theta,\varphi) (\phi_{il}\nabla\phi_{jl}+\nabla\phi_{il}\phi_{jl}) \nabla Y_{lm}

There is an interesting identity:

.. math::

    \int r{\bf \hat r} Y_{lm} \nabla Y_{l'm'} \sin\theta\,\d\theta\,\d\varphi= 0

But it cannot be applied, because we have one more $r$ in the expression.
Nevertheless the term is probably zero, as can be seen when we compare the weak
formulation to the one we got directly from the radial equation.

How Not To Derive The Weak Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


If we forgot that we have to integrate covariantly, this section is devoted
to what happens if we integrate using the coordinate integration. We would get:

.. math::

    \int_0^a {1\over2}R'(x)v'(x)-{1\over r}R'(x)v(x)+ \left(V+{l(l+1)\over2r^2}\right)R(x)v(x)\,\d x = E\int_0^aR(x)v(x)\,\d x

Notice the matrix on the left hand side is not symmetric. There is another way
of writing the weak formulation by applying per-partes to the $R'(r)v(r)$ term:

.. math::

    -\int_0^a{1\over r}R'(x)v(x)\d x=


.. math::

     =\int_0^a{1\over r}R(x)v'(x)\d x -\int_0^a{1\over r^2}R(x)v(x)\d x -\left[{1\over r}R'(x)v'(x)\right]_0^a +\left[{1\over r^2}R'(x)v(x)\right]_0^a

We can use $v(a)=0$ and $R'(a)=0$ to simplify a bit:

.. math::

    -\int_0^a{1\over r}R'(x)v(x)\d x=


.. math::

     =\int_0^a{1\over r}R(x)v'(x)\d x -\int_0^a{1\over r^2}R(x)v(x)\d x +\lim_{r\to0}\left({R'(x)v'(x)\over r}-{R'(x)v(x)\over r^2}\right)

Since $R(x)\sim r^l$ near $r=0$, we can see that for $l\ge3$ the limits
on the right hand side are zero, but for $l=0, 1, 2$ they are not zero and need
to be taken into account. Let's assume $l\ge3$ for now, then our weak formulation looks like:

.. math::

    \int_0^a {1\over2}R'(x)v'(x)+{1\over r}R(x)v'(x)+ \left(V+{l(l+1)\over2r^2}-{1\over r^2}\right)R(x)v(x)\,\d x = E\int_0^aR(x)v(x)\,\d x

or

.. math::

    \int_0^a {1\over2}R'(x)v'(x)+{1\over r}R(x)v'(x)+ \left(V+{(l-2)(l+1)\over2r^2}\right)R(x)v(x)\,\d x = E\int_0^aR(x)v(x)\,\d x

The left hand side is also not symmetric, however we can now take an average of
our both weak formulations to get a symmetric weak formulation:

.. math::

    \int_0^a {1\over2}R'(x)v'(x)+{R(x)v'(x)-R'(x)v(x)\over 2r}+ \left(V+{l(l+1)-1\over2r^2}\right)R(x)v(x)\,\d x =


.. math::

     = E\int_0^aR(x)v(x)\,\d x

Keep in mind, that this symmetric version is only correct for $l\ge3$. For
$l<3$ we need to use our first nonsymmetric version.

As you can see, this is something very different to what we got in the previous
section. First there were lots of technical difficulties and second the final
result is wrong, since it doesn't correspond to the 3D Schroedinger equation.

TODO list
---------


Currently, the code can handle an arbitrary number of equations and solve them
with elements up to the 10th degree. However, the meshes still have to be the
same for every solution component. The code is not $hp$-adaptive yet. These
things will be fixed as time permits.
