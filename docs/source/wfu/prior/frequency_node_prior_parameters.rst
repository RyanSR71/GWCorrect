frequency_node_prior_parameters
===============================

.. code-block:: python

  GWCorrect.wfu.prior.frequency_node_prior_parameters(waveform_generator,asd_data,prior,samples=10000,xi_max=0.318)

Generates gaussian parameters for the xi_0 and delta_xi_tilde truncated gaussian priors by comparing waveforms to detector noise.

.. math::

  \sqrt{\frac{4c^3\xi}{GM}}|\mu(\xi;\vartheta)|-\sqrt{S_n(\xi)}=0

Parameters:
-----------
waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
    bilby waveform generator object
asd_data: numpy.ndarray
    array of amplitude spectral density data; assumes standard LIGO formatting
prior: bilby.core.prior.PriorDict
    bilby prior object to draw BBH parameters from
samples: int, optional, (10000)
    number of waveforms to draw
xi_max: float, optional, (1/pi, 0.318...)
    upper bound on the dimensionless frequency grid

Returns:
--------
mu_down: float
    mean of the xi_0 prior
sigma_down: float
    standard deviation of the xi_0 prior
mu: float
    mean of the delta_xi_tilde prior
sigma: float
    standard deviation of the delta_xi_tilde prior
