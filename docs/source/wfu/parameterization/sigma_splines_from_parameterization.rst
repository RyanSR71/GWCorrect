sigma_splines_from_parameterization
===================================

.. code-block:: python

   GWCorrect.wfu.parameterization.sigma_splines_from_parameterization(parameterization,dimensionless=True,xi_min=0.001,xi_max=0.318,resolution=1000)

Generates the standard deviations of the waveform correction priors and stores them in spline objects from a parameterization file.

.. math::

    \sigma_\mathcal{A}(\xi)=\sqrt{(\overline{\Delta\mathcal{A}_\mu}(\xi))^2+\left(\delta\mathcal{A}_\mu(\xi)\right)^2}

.. math::

    \sigma_\phi(\xi)=\sqrt{(\overline{\Delta\phi_\mu}(\xi))^2+\left(\delta\phi_\mu(\xi)\right)^2}

Parameters:
-----------
parameterization: numpy.ndarray
    parameterization file from GWCorrect.wfu.parameterization.parameterization
dimensionless: bool, optional
    whether or not the output is returned in dimensionless frequency units
    Default: True
xi_min: float, optional
    if dimensionless is True, this is the lower bound on the dimensionless frequency grid
    default: 0.001
xi_max: float, optional
    if dimensionless is True, this is the upper bound on the dimensionless frequency grid
    default: 1/pi, 0.318...
resolution: int, optional
    if dimensionless is True, this is the number of points in the dimensionless frequency grid
    default: 1000

Returns:
--------
sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dA prior
sigma_dphi_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dphi prior
