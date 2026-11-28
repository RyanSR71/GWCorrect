ringdown_frequency
==================
.. code-block:: python

  GWCorrect.ppE.utils.ringdown_frequency(strain,frequency_array, total_mass)

Numerically finds the ringdown frequency, f_RD, given a frequency domain waveform. The ringdown frequency is as defined in: https://arxiv.org/pdf/1508.07253

.. math::

  \left(\frac{\partial^2\phi}{\partial f^2}\right)_{f=f_\mathrm{RD}}=0

Parameters:
-----------
fd_waveform : FrequencySeries
  The frequency domain waveform of the coalescence.
total_mass : float
  The total mass of the binary, in solar masses.

Returns:
--------
f_RD : float
    The ringdown frequency, in Hz.
