double_gaussian
===============
.. code-block:: python

  GWCorrect.wfu.utils.double_gaussian(x, N1, mu1, sigma1, N2, mu2, sigma2)

Double Gaussian function. Used with GWCorrect.wfu.utils.gaussian_parameters_from_A_ASD_solutions

.. math::

  N_1e^{-\frac{1}{2}\left(\frac{x-\mu_1}{\sigma_1}\right)^2}+N_2e^{-\frac{1}{2}\left(\frac{x-\mu_2}{\sigma_2}\right)^2}

Parameters:
-----------
x: numpy.ndarray
    input array
N1: float
    first scaling constant
mu1: float
    mean of the first Gaussian
sigma1: float
    standard deviation of the first Gaussian
N2: float
    second scaling constant
mu2: float
    mean of the second Gaussian
sigma2: float
    standard deviation of the second Gaussian

Returns:
--------
output_array: numpy.ndarray
    output double Gaussian corresponding to the input array and parameters
