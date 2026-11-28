apply_tapering_to_amplitude
===========================
.. code-block:: python

  GWCorrect.ppE.utils.apply_tapering_to_amplitude(amplitude,e)

Applies the one-sided Planck taper window to a frequency domain amplitude. 

Parameters:
-----------
amplitude: numpy.ndarray
  frequency domain amplitude
e: float
  taper parameter
                           
Returns:
--------
tapered_amplitude: numpy.ndarray
  input amplitude with taper window applied
