keplerian_velocity
==================
.. code-block:: python

  GWCorrect.ppE.utils.keplerian_velocity(freqs, total_mass, ppe_b)

Computes the Keplerian velocty. This function is defined assuming geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses. Likewise for if the total mass is given seconds.

.. math::

  u\equiv(\pi Mf)^{b/3}

Parameters:
-----------
freqs: numpy.ndarray
  array of frequencies in geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses, likewise for if the total mass is given seconds
total_mass: float
  total binary mass in geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses, likewise for if the total mass is given seconds
ppe_b: float
  ppE PN order parameter

Returns:
--------
u: int
  Keplerian velocity
