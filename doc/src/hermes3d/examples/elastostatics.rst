Elastostatics
=============

**Git reference:** Example `elastostatics <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes3d/examples/elastostatics>`_.

This example deals with equations of linear elasticity inside an L-shaped domain. Elastostatics studies 
linear elastic deformations under the conditions of equilibrium where all forces acting on the elastic 
body sum to zero, and  displacements are not a function of time. 

.. index::
   single: mesh; fixed
   single: problem; elliptic, linear, symmetric

The governing equations have the form:

.. math::
   :nowrap:
   :label: elastostatics

   \begin{eqnarray*}
   \sigma_{ji,j} + F_i & = & 0 \hbox{ in }\Omega, \\ \nonumber
   \epsilon_{ij}       & = & \frac{1}{2}(u_{j,i} + u_{i,j}),   \\
   \sigma_{i,j}        & = & C_{ijkl} \, \epsilon_{kl}.
   \end{eqnarray*}

Here the subscript $\cdot_{,j}$ indicates $\partial{\cdot}/\partial x_j$, $\sigma_{ji,j}$ is the 
stress tensor, $\epsilon_{ij}$ is the strain (deformation), $u_i$ is the displacement,
$C_{ijkl}$ is the forth-order stiffness tensor. By Einstein summation convention, 
the $3^{rd}$ equation of :eq:`elastostatics` represent the following: 

.. math::
   :nowrap:
   :label: elasto-sum

   \begin{eqnarray*}
   C_{ijkl} \, \epsilon_{kl} & = & \sum_{k,l=1}^3 C_{ijkl} \, \epsilon_{kl},
   \end{eqnarray*}

where $1 \le i, j, k, l \le 3$.

.. image:: elastostatics/elasto-statics-domain.png

The domain of interest is an L-shaped beam equipped with 
zero Dirichlet boundary conditions: $u_1 = u_2 = u_3 = 0$ on all five boundary faces (${\Gamma}_u$) 
except the left-most vertical one (${\Gamma}_F$), where an external force $F$ is applied::

        // Boundary condition types.
        BCType bc_types_x(int marker)
        {
          return BC_NATURAL;
        }

        BCType bc_types_y(int marker)
        {
          return BC_NATURAL;
        }

        BCType bc_types_z(int marker)
        {
          return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL;
        }


The stiffness tensor $C_{ijkl}$ is constant and symmetric:

.. math::
   :nowrap:
   :label: elasto-stress

   \begin{eqnarray*}
   \sigma_{ij} & = & \lambda \delta_{ij} \epsilon_{kk} + 2\mu\epsilon_{ij}, \\ \nonumber
   \lambda     & = & \frac{E\nu}{(1+\nu)(1-2\nu)},                          \\
   \mu         & = & \frac{E}{2(1+\nu)}. 
   \end{eqnarray*}

Here $\lambda, \mu$ are the Lame constants, $E$ is the Young modulus, $\nu$ is the Poisson ratio. 
In our example, $E = 200 \times 10^9$ Gpa and $\nu = 0.3.$ 

Substituting :eq:`elasto-stress` back into :eq:`elastostatics`, we obtain:
 
.. math::
   :nowrap:
   :label: elasto-navier

   \begin{eqnarray*}
   \mu u_{i,jj}  + (\mu + \lambda)u_{j,ij} + F_i & = & 0,              \\ \nonumber
   \hbox{ or }           & \, & \\                                      
   \mu \Delta{u} + (\mu + \lambda) \mathsf{grad} \, \mathsf{div} u  + F & = & 0.
   \end{eqnarray*}

The corresponding weak formulation is as follows:

.. math::
   :nowrap:
   :label: elasto-statics-form

   \begin{eqnarray*}
   \int_{\Omega} (\lambda + 2\mu) u_{i} \, v_{i} + \mu u_{j} \, v_{j} + \mu u_{k} \, v_{k} \quad 
   +\quad \int_{\Omega} \lambda u_{i} \,  v_{j} + \mu u_{j} \, v_{i} \quad
   +\quad \int_{\Omega} \lambda u_{i} \,  v_{k} + \mu u_{k} \, v_{i}
     &  = & 0, \\ \nonumber
   \int_{\Omega} \mu u_{i} \, v_{i} + (\lambda + 2\mu) u_{j} \, v_{j} + \mu u_{k} \, v_{k} \quad
   +\quad \int_{\Omega} \lambda u_{j} \,  v_{k} + \mu u_{k} \, v_{j}
     &  = & 0, \\
   \int_{\Omega} \mu u_{i} \, v_{i} + \mu u_{j} \, v_{j} + (\lambda + 2\mu) u_{k} \, v_{k} 
     &  = & \int_{\Gamma_F} F_i v. \nonumber
   \end{eqnarray*}

Here is the code for the weak forms::

    template<typename real, typename scalar>
    scalar bilinear_form_0_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dx_b_dy_c_dz<real, scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
    }
      
    template<typename real, typename scalar>
    scalar bilinear_form_0_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dx_b_dy_c_dz<real, scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
    }

    template<typename real, typename scalar>
    scalar bilinear_form_0_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dudx_dvdy_b_dudy_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
    }

    template<typename real, typename scalar>
    scalar bilinear_form_0_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dudx_dvdz_b_dudz_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
    }

    template<typename real, typename scalar>
    scalar surf_linear_form_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return 0.0;
    }

    template<typename real, typename scalar>
    scalar bilinear_form_1_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dx_b_dy_c_dz<real, scalar>(mu, lambda + 2*mu, mu, n, wt, u, v, e);
    }

    template<typename real, typename scalar>
    scalar bilinear_form_1_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dudy_dvdz_b_dudz_dvdy<real, scalar>(lambda, mu, n, wt, v, u, e);
    }

    template<typename real, typename scalar>
    scalar surf_linear_form_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return 0.0;
    }

    template<typename real, typename scalar>
    scalar bilinear_form_2_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      return int_a_dx_b_dy_c_dz<real, scalar>(mu, mu, lambda + 2*mu, n, wt, u, v, e);
    }

    template<typename real, typename scalar>
    scalar surf_linear_form_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
    {
      scalar res = 0.0;
      for (int i = 0; i < n; i++) res += wt[i] * (f * v->fn[i]);
      return res;
    }

Solution graph:

.. image:: elastostatics/elasto-statics-sln.png

