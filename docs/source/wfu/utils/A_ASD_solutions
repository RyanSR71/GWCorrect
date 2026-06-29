A_ASD_solutions
===============

.. code-block:: python

  GWCorrect.wfu.utils.A_ASD_solutions(waveform_generator,asd_data,prior,samples,desc)

Computes the distributions in dimensionless frequency where waveforms cross detector noise.

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
samples: int
    number of waveforms to draw
desc: string
    description of what this function is doing; this will be printed by the logger while this function runs

Returns:
--------
lower_xis: list
    list of dimensionless frequency values where a frequency domain waveform enters the detectable region of a detector
upper_xis: list
    list of dimensionless frequency values where a frequency domain waveform leaves the detectable region of a detector
