phase_difference_from_parameters
================================
.. code-block:: python

  GWCorrect.wfu.utils.phase_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,
                                                       phi_tildes,sigma_phi_spline,mass_1,mass_2,gamma)

Generates the amplitude correction from the frequency nodes and rescaled amplitude correction parameters

.. math::

  \Delta\phi_\mathrm{SI}=\frac{1}{2\gamma N}\int_{f-\gamma N}^{f+\gamma N}\Delta\phi_\mathrm{1D}(f^\prime;\{f_k,\tilde\varphi_k\sigma_\phi(\xi_k)\})df^\prime

Parameters:
-----------
frequency_array: numpy.ndarray
    array of SI frequency points over which to generate the phase correction
dimensionless_frequency_nodes: numpy.ndarray
    array of frequency nodes in dimensionless frequency nodes
alpha_tildes: list
    list of rescaled amplitude correction parameters
sigma_phi_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the total standard deviation of the phase correction prior
mass_1: float
    primary mass in solar masses
mass_2: float
    secondary mass in solar masses
gamma: float
    smoothing parameter for the smooth interpolation function

Returns:
--------
phase_difference: numpy.ndarray
    array of phase differences with the same shape as the input frequency_array
