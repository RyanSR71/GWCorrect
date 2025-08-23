match
=====

.. code-block:: python

  GWCorrect.ppE.prior.match(waveform_1,waveform_2,injection)

Computes the normalized match between two waveforms.

.. math::

  \mathrm{match}(\tilde{h}_1,\tilde{h}_2)=\underset{\Delta t_c,\Delta\phi_c}{\mathrm{max}}\Re\left[\frac{\int\hat{h}_1^*(f)\hat{h}_2(f)\exp(-2\pi if\Delta t_c+i\Delta\phi_c)df}{\sqrt{\int \hat{h}_1^*(f)\hat{h}_1(f)df}\sqrt{\int \hat{h}_2^*(f)\hat{h}_2(f)df}}\right]

Parameters:
-----------
waveform_1: bilby.gw.WaveformGenerator
    bilby waveform generator object
waveform_2: bilby.gw.WaveformGenerator
    bilby waveform generator object
injection: dict
    dictionary of source parameters to be put into the waveform generators 

Returns:
--------
match: float
    the match between the two waveforms; ranges from 0 to 1
