

The linear chain
=================

Let us consider the Hamiltonian of a periodic linea chain:

.. math::

    H = t \sum_{i} \ket{i}\bra{i+1} + {\rm h.c}.

The eigenvalues of this system :math:`\varepsilon_n =-2t \cos(k_na)` and due to the periodic boundary conditions the momentum is constrain to the 
following values  :math:`k_n=\frac{2\pi n}{aN}`. If we impose periodic boundary conditions, 
the 

The density of states
---------------------

The density of states (DOS) at a given energy :math:`\rho(\varepsilon)` is defined as the number of states within the
interval :math:`(\varepsilon,\varepsilon+d\varepsilon)`. In the periodic linear chain, the density of states could be computed by
counting the allowed states within a interval, which can be achieved through the following function

.. math::

    \rho(\varepsilon) = \frac{1}{N} \int_{\rm all energies}  {\rm Tr} [\delta(H-E)] {\rm d}E.

The  expression above is general and relies on the property of the delta. Alternatively, for periodic systems, one can also integrate the density of momentums,
since it is univocally linked to the density of states. To make this simple, let us consider the limit of infinite sytem, where we can approximate :math:`x=n/N` as a real number 
within the :math:`[0,1)` interval and :math:`ka= 2\pi x` is also a continuum value.  Using the dispersion relation, we can compute the density of energy per momentum 

.. math::

    \frac{d\varepsilon_k}{dk}  =2t a  \sin(ka)  

and through proper inversion and by exploiting the relation between momentum and state index, also the ammount of states per energy which is nothing but the DOS

.. math::

    \rho(\varepsilon) &= \frac{dx(\varepsilon) }{d\varepsilon}, \\
                      &= \frac{1}{2\pi \sqrt{4t^2- \varepsilon^2} }.
    
The kernel polynomial method 
-----------------------------

If we choose the hopping energy in the periodic linear change to be :math:`t=1/2` then the spectrum is bounded within the :math:`(-1,1)` interval and this make it suitable for KPM expansion.
As described in the kpm_module computing the density of states in the KPM approach involves calculating the moments

.. math::

    \mu_m  = \frac{1}{N} Tr[ T_m(H) ]

where :math:`T_m(x)\equiv = \cos(m \arccos(x))`. If we choose the eigenvector basis, the above is simply

.. math::

    \mu_m  =\frac{1}{N} \sum_{n} \cos\left[m \arccos\left(\cos\frac{2\pi n}{N}\right)\right] 

In the limit of infinite system, the sum can be replace by an integral 

.. math::

    \mu_m = \int_{0}^1 dx \cos\left(2\pi  m x\right)= \delta_{m0}

which after replacing in the chebyshev sum leads to

.. math::

    \rho(\varepsilon) = \frac{1}{\pi \sqrt{1- x^2}

fin





