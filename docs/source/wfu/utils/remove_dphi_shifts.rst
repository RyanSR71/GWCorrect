remove_dphi_shifts
==================
.. code-block:: python

  GWCorrect.wfu.utils.remove_dphi_shifts(dimensionless_frequency_grid,injection,n,sigma_dphi_spline,
                                         asd_data,waveform_generator,xi_max=0.318,gamma=0.025)

Takes an injection of binary black hole and phase correction parameters and removes any time domain shifts from the phase deviation.

.. math::

  \Delta\phi_\mathrm{SI}^{\mathrm{aligned}}(\xi;\Phi)=\Delta\phi_\mathrm{SI}(\xi;\Phi)-\frac{2\pi c^3\xi t_0}{GM}-\phi_0

Parameters:
-----------
dimensionless_frequency_grid: numpy.ndarray
    array of dimensionless frequencies to evaluate the phase difference over
injection: dict
    dictionary of binary black hole source parameters and phase correction parameters
n: int
    number of frequency nodes excluding 0
sigma_dphi_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dphi prior
asd_data: numpy.ndarray
    array of amplitude spectral density data; assumes this is in the standard LIGO format
waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
    bilby waveform generator object
xi_max: float, optional, (1/pi, 0.318...)
    upper bound on the dimensionless frequency band
gamma: float, optional, (0.025)
    smoothing parameter
  
Returns:
--------
phase_difference_no_shifts: numpy.ndarray
    input phase difference from injection with overall time and phase shifts removed; same shape as input dimensionless frequency grid
