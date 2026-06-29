delta_t_01
====================

.. code-block:: python

  GWCorrect.wfu.utils.delta_t_01(xi_0,delta_xi_tilde,total_mass,n,xi_max)

delta_t_01 conversion function.

.. math::

  \Delta t_{01} = \frac{nGM}{c^3(\xi_\mathrm{max}-\xi_0)\delta\tilde\xi}

Parameters:
-----------
xi_0: float
    starting frequency node parameter
delta_xi_tilde: float
    frequency node spacing parameter
total_mass: float
    total system mass (m_1+m_2)
n: int
    number of frequency nodes excluding 0
xi_max: float
    upper bound on the dimensionless frequency band

Returns:
--------
delta_t_01: float
    delta_t_01 value
