dh_SI
=====
.. code-block:: python

  GWCorrect.wfu.utils.dh_SI(dimensionless_frequency_grid,dimensionless_frequency_nodes,dhs,sigma_dh_spline,gamma=0.025)

Smooth interpolation function for the waveform correction parameters.

Parameters:
-----------
dimensionless_frequency_grid: numpy.ndarray
    array of mass-weighted dimensionless frequencies
dimensionless_frequency_nodes: list or numpy.ndarray
    list of dimensionless frequency nodes
dhs: list or numpy.ndarray
    rescaled waveform correction parmaters, either dA or dphi
sigma_dh_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dh prior, either dA or dphi
gamma: float, optional, (0.025)
    smoothing parameters

Returns:
--------
dh_smooth_interpolation: numpy.ndarray
    dA or dphi smooth interpolation array corresponding to the input dimensionless frequency grid
