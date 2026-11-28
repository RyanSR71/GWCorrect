outer_frequency_node_prior_parameters
=====================================

.. code-block:: python

  GWCorrect.wfu.prior.outer_frequency_node_prior_parameters(waveform_generator,asd_data,prior,
  samples=10000,xi_max=0.318)

Generates gaussian parameters for the xi_0 and xi_n truncated gaussian priors by comparing waveforms to detector noise.

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
minimum_frequency: float
    lower bound on the SI frequency band
maximum_frequency: float
    upper bound on the SI frequency band
samples: int, optional, (10000)
    number of waveforms to draw

Returns:
--------
mu_0: float
    mean of the lower_xis distribution
sigma_0: float
    standard deviation of the lower_xis distribution
mu_n: float
    mean of the upper_xis distribution
sigma_n: float
    standard deviation of the upper_xis distribution
