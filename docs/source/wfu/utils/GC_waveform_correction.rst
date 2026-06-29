GC_waveform_correction
======================

.. code-block:: python

  GWCorrect.wfu.utils.BBH_waveform_correction(frequency_array,xi_0,delta_xi_tilde,dAs,dphis,sigma_dA_spline,
                                              sigma_dphi_spline,mass_1,mass_2,xi_max,gamma)

Computes (1+dA)exp(idphi) from the frequency node and correction parameters.

.. math::

  (1+\Delta\mathcal{A}_\mathrm{SI}(f;\mathrm{A}))\exp(i\Delta\phi_\mathrm{SI}(f;\Phi))

Parameters:
-----------
frequency_array: numpy.ndarray
    frequency grid over which to generate the waveform correction
xi_0: float
    starting frequency node parameter
delta_xi_tilde: float
    frequency node spacing parameter
dAs: list or numpy.ndarray
    list of dA (alpha_tilde) parameters
dphis: list or numpy.ndarray
    list of dphi (varphi_tilde) parameters
sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dA prior
sigma_dphi_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dphi prior
mass_1: float
    mass of the primary body in the binary
mass_2: float
    mass of the secondary body in the binary
xi_max: float
    upper bound on the dimensionless frequency grid
gamma: float
    smoothing parameter for the smooth interpolation function

Returns:
--------
waveform_correction: numpy.ndarray
    (1+dA)exp(idphi) array to multiply to frequency domain gravitational-wave strain; same shape as frequency_array input
