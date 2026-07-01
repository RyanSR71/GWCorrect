conversion
==========

.. code-block:: python

  GWCorrect.wfu.prior.conversion(parameters,n=None,xi_max=0.318,total_mass=None,sigma_dA_spline=None)

Conversion function needed to use the total_mass, delta_t_01, and epsilon_alpha constraint priors. Used with the bilby.core.prior.PriorDict() function

.. math::

  \Delta t_{01}\equiv\frac{nGM}{c^3(\xi_\mathrm{max}-\xi_0)\delta\tilde\xi}

.. math::

  \epsilon_\alpha = n-\sum_{k=1}^n\theta(\tilde\alpha_k\sigma_\mathcal{A}(\xi_k)+1)

Parameters:
------------
input_parameters: dict
    dictionary of binary black hole source parameters and waveform correction parameters
n: int, optional, (None)
        number of frequency nodes excluding 0; if given, delta_t_01 and epsilon_alpha will be calculated
xi_max: float, optional, (1/pi, 0.318...)
    upper bound on the dimensionless frequency band
total_mass: float, optional, (None)
    if mass parameters are not being sampled in the prior, set the fixed total mass value here
sigma_dA_spline: scipy.interpolate._cubic.CubicSpline, optional, (None)
    scipy cubic spline object encoding the standard deviation of the dA prior; if given, epsilon_alpha will be calculated

Returns:
--------
parameters: dict
    input parameters plus total mass, delta_t_01, and epsilon_alpha (assuming all necessary parameters are present)
