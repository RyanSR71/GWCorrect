convert_to_SI_frequency
==================================

.. code-block:: python

  GWCorrect.wfu.conversion.convert_to_SI_frequency(dimesionless_frequencies,total_mass)

Converts mass-weighted dimensionless frequency to SI frequency.
  
.. math::

  f=\frac{c^3\xi}{GM}

Parameters:
-----------
dimensionless_frequencies: float or numpy.ndarray
    mass-weighted dimensionless frequency values
total_mass: float
    total binary mass to weight the frequencies

Returns:
--------
SI_frequency: float or numpy.ndarray
    SI frequencies of the same shape as the input dimensionless frequencies
