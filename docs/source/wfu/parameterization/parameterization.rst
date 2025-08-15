parameterization
================

.. code-block:: python

   GWCorrect.wfu.parameterization.parameterization(hf1,hf2,prior,nsamples,
   spline_resolution=0.5,npoints=None,psd_data=None,correction_parameter=1e-4,
   ref_amplitude=None,polarization='plus')

Generates samples of waveform differences between two approximants and parameterizes the data.

Parameters:
-----------
hf1: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object
hf2: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object
prior: bilby.core.prior.dict.PriorDict
   bilby prior object
nsamples: int
   number of draws of waveform uncertainty desired
spline_resolution: float, optional, (0.5)
   fraction of the frequency grid length to use as the number of spline nodes
npoints: int, optional, (None)
   length of the desired frequency grid; if None, this will be set automatically by the frequency grid
psd_data: numpy.ndarray, optional, (None)
   array containing the power spectral density data and their corresponding frequencies; if None, psd_data will be set to `GW170817 PSD data <https://dcc.ligo.org/ligo-p1900011/public>`_ 
correction_parameter: float, optional, (0.0001)
   fraction of maximum amplitude to cut off the amplitude at
ref_amplitude: numpy.ndarray, optional, (None)
   reference amplitude for residual phase calculation
polarization: string, optional, ('plus')
   polarization of the strain data (plus or cross)
  
Returns:
--------
parameterized_data: numpy.ndarray
   table containing the following:
      frequency_grid: numpy.ndarray
         frequencies corresponding to the frequency parameters specified
      frequency_nodes: numpy.ndarray
         frequency nodes for the splines
      dA_parameters: numpy.ndarray
         amplitude difference spline parameters
      dphi_parameters: numpy.ndarray
         phase difference spline parameters
      injection: dictionary
         source parameters injected into the waveform generators
