Mathematical Background
-----------------------

Main strengths of Hermes are 

 * Robust PDE-independent adaptivity algorithms. 
 * Adaptivity algorithms for time-dependent problems.
 * Monolithic discretization of multiphysics problems.

The following list describes the above in more detail:

* **Mature hp-adaptivity algorithms**. Hermes puts a major emphasis on error control and automatic adaptivity. Practitioners know well how painful it is to use automatic adaptivity in conjunction with standard lower-order approximations such as linear or quadratic elements - the error decreases somehow during a few initial adaptivity steps, but then it slows down and it does not help to invest more unknowns or CPU time. This is typical for low-order methods. In contrast to this, the exponentially-convergent adaptive *hp*-FEM and *hp*-DG do not have this problem - the error drops steadily and fast during adaptivity all the way to the desired accuracy. The following graph shows typical convergence rates of *h*-FEM with linear elements, *h*-FEM with quadratic elements, and *hp*-FEM on a log-log scale:

.. image:: hermes2d/img/intro/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: Typical convergence curves of FEM with linear and quadratic elements and hp-FEM

Same graphs as above but now in terms of CPU time:

.. image:: hermes2d/img/intro/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

* **Wide applicability**. Hermes is PDE-independent. Standard FEM codes are designed to solve some narrow class(es) of PDE problems (such as elliptic equations, fluid dynamics, electromagnetics etc.). Hermes does not employ any technique or algorithm that would limit its applicability to some particular class(es) of PDE problems. Automatic adaptivity is guided by a universal computational a-posteriori error estimate that works in the same way for any PDE. Of course this does not mean that the algorithms perform equally well on all PDE - some equations simply are more difficult to solve than others. However, Hermes allows you to tackle an arbitrary PDE or multiphysics PDE system and add your own equation-specific extensions if necessary. Visit the `hp-FEM group home page <http://hpfem.org/>`_ and especially the `gallery <http://hpfem.org/gallery/>`_ to see numerous examples.

.. image:: hermes2d/img/intro/ns.jpg
   :align: center
   :width: 650
   :height: 300
   :alt: Image of incompressible viscous flow.


* **Arbitrary-level hanging nodes**. Hermes has a unique original methodology for handling irregular meshes with arbitrary-level hanging nodes. This means that extremely small elements can be adjacent to very large ones. When an element is refined, its neighbors are never split forcefully as in conventional adaptivity algorithms. It is well known that approximations with one-level hanging nodes are more efficient compared to regular meshes. However, the technique of arbitrary-level hanging nodes brings this to a perfection.

.. image:: hermes2d/img/intro/ord_2d_c.png
   :align: center
   :width: 370
   :height: 350
   :alt: Illustration of arbitrary-level hanging nodes.

.. ######
    .. image:: hermes2d/img/intro/mixer-mesh.png
       :align: right
       :width: 300
       :height: 300
       :alt: Illustration of arbitrary-level hanging nodes.

    .. raw:: html

       <hr style="clear: both; visibility: hidden;">

* **Multimesh hp-FEM**. Various physical fields or solution components in multiphysics problems can be approximated on individual meshes, combining quality $H^1$, $H(curl)$, $H(div)$, and $L^2$ conforming higher-order elements. Due to a unique original methodology, no error is caused by operator splitting, transferring data between different meshes, and the like. The following figure illustrates a coupled problem of heat and moisture transfer in massive concrete walls of a nuclear reactor vessel. 

.. image:: hermes2d/img/intro/hm-sln-frame.png
   :align: left
   :width: 480
   :alt: Illustration of multimesh hp-FEM.

.. image:: hermes2d/img/intro/hm-mesh-frame.png
   :align: right
   :width: 480
   :alt: Illustration of multimesh hp-FEM.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

* **Dynamical meshes for time-dependent problems**. In time-dependent problems, different physical fields or solution components can be approximated on individual meshes that evolve in time independently of each other. Due to a unique original methodology, no error is caused by transfering solution data between different meshes and time levels. No such transfer takes place in the multimesh *hp*-FEM - the discretization of the time-dependent PDE system is monolithic. 

.. image:: hermes2d/img/intro/flame.jpg
   :align: center
   :width: 700
   :height: 360
   :alt: Adaptive hp-FEM with dynamical meshes for a flame propagation problem. 
