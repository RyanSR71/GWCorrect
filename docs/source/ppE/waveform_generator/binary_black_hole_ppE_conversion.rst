binary_black_hole_ppE_conversion
================================

.. code-block:: python

   GWCorrect.ppE.waveform_generator.binary_black_hole_ppE_conversion(parameters)

EXPERIMENTAL: Converts set of parameters into the parameters needed by MultippECorrectionModel. 

Parameters:
---------------
parameters: dict
    set of binary black hole parameters and correction parameters

Returns:
--------
converted_parameters: dict
    input parameters, but with added parameters needed by the source model
added_keys: list
    list of parameter names that were added to the output parameter dictionary
