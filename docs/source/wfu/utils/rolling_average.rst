rolling_average
===============
.. code-block:: python

  GWCorrect.wfu.utils.rolling_average(data,gamma)

Rolling average function.

.. math::

  \int_{x-\gamma N}^{x+\gamma N} f(x^\prime)dx^\prime

Parameters:
-----------
data: numpy.ndarray
    input array to be averaged
gamma: float
    smoothing parameter

Returns:
--------
new_data: numpy.ndarray
    averaged input data
