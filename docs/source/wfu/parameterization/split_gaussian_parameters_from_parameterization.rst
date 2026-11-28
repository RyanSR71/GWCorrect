split_gaussian_parameters_from_parameterization
===============================================

.. code-block:: python

   GWCorrect.wfu.parameterization.split_gaussian_parameters_from_parameterization(parameterization,
   dimensionless=True,xi_min=0.001,xi_max=0.318,resolution=1000)

Takes all of the sets in a parameterized waveform difference matrix and takes the mean and standard deviation of amplitude and phase difference.

Parameters:
-----------
parameterization: numpy.ndarry
   parameterization matrix
dimensionless: bool, optional, (True)
   whether or not the output is returned in dimensionless frequency units
xi_min: float, optional, (0.001)
   if dimensionless is True, this is the lower bound on the dimensionless frequency grid
xi_max: float, optional, (1/pi, 0.318...)
   if dimensionless is True, this is the upper bound on the dimensionless frequency grid
resolution: int, optional, (1000)
   if dimensionless is True, this is the number of points in the dimensionless frequency grid

Returns:
--------
dimensionless_frequency_grid: numpy.ndarray
   array of dimensionless frequency values corresponding to the output arrays
dA_minus: numpy.ndarray
   standard deviation below the median amplitude difference
median_dA: numpy.ndarray
   median amplitude difference values
dA_plus: numpy.ndarray
   standard deviation above the median amplitude difference
dphi_minus: numpy.ndarray
   standard deviation below the median phase difference
median_dphi: numpy.ndarray
   median phase difference values
dphi_plus: numpy.ndarray
   standard deviation above the median phase difference    
