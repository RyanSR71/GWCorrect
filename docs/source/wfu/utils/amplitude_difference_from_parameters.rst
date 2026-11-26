amplitude_difference_from_parameters
====================================
.. code-block:: python

  GWCorrect.wfu.utils.amplitude_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,
                                                           alpha_tildes,sigma_A_spline,mass_1,mass_2,gamma)

Generates the amplitude correction from the frequency nodes and rescaled amplitude correction parameters

.. math::

  \Delta\mathcal{A}_\mathrm{SI}=\frac{1}{2\gamma N}\int_{f-\gamma N}^{f+\gamma N}\Delta\mathcal{A}_\mathrm{1D}(f^\prime;\{f_k,\tilde\alpha_k\sigma_\mathcal{A}(\xi_k)\})df^\prime

Parameters:
-----------
frequency_array: numpy.ndarray
    array of SI frequency points over which to generate the amplitude correction
dimensionless_frequency_nodes: numpy.ndarray
    array of frequency nodes in dimensionless frequency nodes
alpha_tildes: list
    list of rescaled amplitude correction parameters
sigma_A_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the total standard deviation of the amplitude correction prior
mass_1: float
    primary mass in solar masses
mass_2: float
    secondary mass in solar masses
gamma: float
    smoothing parameter for the smooth interpolation function

Returns:
--------
amplitude_difference: numpy.ndarray
    array of amplitude differences with the same shape as the input frequency_array
