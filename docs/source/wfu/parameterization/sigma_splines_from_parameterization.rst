sigma_splines_from_parameterization
===================================

.. code-block:: python

   GWCorrect.wfu.parameterization.sigma_splines_from_parameterization(parameterization,
   dimensionless=True,xi_min=0.001,xi_max=0.318,resolution=1000)

Generates the standard deviations of the waveform correction priors and stores them in spline objects from a parameterization file.

.. math::

    \overline{\delta\mathcal{A}_\mu}=\sqrt{\left(1-2/\pi\right)\left(\delta\mathcal{A}_\mu^+-\delta\mathcal{A}_\mu^-\right)^2+\delta\mathcal{A}_\mu^+\delta\mathcal{A}_\mu^-}

.. math::

    \overline{\delta\phi_\mu}=\sqrt{\left(1-2/\pi\right)\left(\delta\phi_\mu^+-\delta\phi_\mu^-\right)^2+\delta\phi_\mu^+\delta\phi_\mu^-}

.. math::

   \sigma_\mathcal{A}(\xi)=\sqrt{\left(\overline{\Delta\mathcal{A}_\mu}(\xi)\right)^2+\left(\overline{\delta\mathcal{A}_\mu}(\xi)\right)^2}

.. math::

   \sigma_\phi(\xi)=\sqrt{\left(\overline{\Delta\phi_\mu}(\xi)\right)^2+\left(\overline{\delta\phi_\mu}(\xi)\right)^2}

Parameters:
-----------
parameterization: numpy.ndarray
    parameterization file from GWCorrect.wfu.parameterization.parameterization
dimensionless: bool, optional, (true)
    whether or not the output is returned in dimensionless frequency units
xi_min: float, optional, (0.001)
    if dimensionless is True, this is the lower bound on the dimensionless frequency grid
xi_max: float, optional, (1/pi, 0.318...)
    if dimensionless is True, this is the upper bound on the dimensionless frequency grid
resolution: int, optional, (1000)
    if dimensionless is True, this is the number of points in the dimensionless frequency grid

Returns:
--------
sigma_A_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the amplitude correction prior
sigma_phi_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the phase correction prior
