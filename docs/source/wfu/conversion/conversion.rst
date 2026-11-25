conversion
==========

.. code-block:: python

  GWCorrect.wfu.conversion.conversion(parameters,n=None,fixed_total_mass=None,xi_max=None,
                                      sigma_A_spline=None,gamma=0.025,duration=None)

Conversion function to generate constrained parameters: total_mass, epsilon_f, and epsilon_alpha

.. math::

  F_{k\ell}=\theta\left(\frac{4(1-\delta_{k0}\delta_{\ell n}(1-n))}{T_\mathrm{obs}}-\frac{c^3(\xi_\ell-\xi_k)}{GM}\right)

.. math::

  \epsilon_{f}\equiv(\delta_{k0}\delta_{\ell n}+\delta_{k,\ell-1})F^{k\ell}

.. math::

  \epsilon_\alpha=\sum_{k=1}^n\theta(-\alpha_k-1)

Parameters:
-----------
input_parameters: dict
    dictionary of BBH source and correction parameters
n: int, optional, (None)
    number of frequency nodes excluding zero; required for epsilon_f and epsilon_alpha calculations
    default: None
fixed_total_mass: float, optional, (None)
    total_mass value for epsilon_f and epsilon_alpha calculations if the total mass is not being sampled in the prior
sigma_A_spline: scipy.interpolate._cubic.CubicSpline, optional, (None)
    scipy cubic spline object encoding the standard deviation of the amplitude correction prior; required for epsilon_alpha calculation
gamma: float, optional, (0.025)
    smoothing parameter for the smooth interpolation function
duration: float, optional, (None)
    observation time duration of the gravitational-wave signal; required for epsilon_f calculation

Returns:
--------
parameters: dict
    dictionary of input parameters plus total mass, epsilon_f, and/or epsilon_alpha; this depends on kwargs
