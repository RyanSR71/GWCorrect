gaussian
========
.. code-block:: python

  GWCorrect.wfu.utils.gaussian(x,N,mu,sigma)

Gaussian function. Used with GWCorrect.wfu.utils.gaussian_parameters_from_A_ASD_solutions

.. math::

  Ne^{-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2}

Parameters:
-----------
x: numpy.ndarray
    input array
N: float
    scaling constant
mu: float
    mean of the Gaussian
sigma: float
    standard deviation of the Gaussian

Returns:
--------
output_array: numpy.ndarray
    output Gaussian corresponding to the input array
