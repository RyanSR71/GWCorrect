conversion
==========

.. code-block:: python

  GWCorrect.wfu.prior.conversion(parameters,n=None,xi_max=0.318)

Conversion function needed to use the total_mass constraint prior and the delta_t_01 constraint prior (optional). Used with the bilby.core.prior.PriorDict() function

.. math::

  \Delta t_{01}\equiv\frac{nGM}{c^3(\xi_\mathrm{max}-\xi_0)\delta\tilde\xi}

 Parameters:
------------
input_parameters: dict
    dictionary of binary black hole source parameters and waveform correction parameters
n: int, optional, (None)
    number of frequency nodes excluding 0; if given, delta_t_01 will be calculated
xi_max: float, optional, (1/pi, 0.318)
    upper bound on the dimensionless frequency band

Returns:
--------
parameters: dict
    input parameters plus total mass and delta_t_01 (if n is not None)
