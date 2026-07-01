convert_to_dimensionless_frequency
==================================

.. code-block:: python

  GWCorrect.wfu.conversion.convert_to_dimensionless_frequency(frequencies,total_mass)

Converts SI frequency to mass-weighted dimensionless frequency.

.. math::

  \xi=\frac{GMf}{c^3}

Parameters:
-----------
frequencies: float or numpy.ndarray
    SI frequency values
total_mass: float
    total binary mass to weight the frequencies

Returns:
--------
dimensionless_frequency: float or numpy.ndarray
    dimensionless frequencies of the same shape as the input frequencies
