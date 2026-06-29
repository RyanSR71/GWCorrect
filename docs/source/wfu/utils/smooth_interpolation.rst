smooth_interpolation
====================

.. code-block:: python

  GWCorrect.wfu.utils.smooth_interpolation(full_grid,nodes,parameters,gamma)

Smooth interpolation function.

.. math::

  f_\mathrm{SI}(x;\lambda)=\frac{1}{2\gamma N}\int_{x-\gamma N}^{x+\gamma N}f_\mathrm{1D}(x^\prime;\{x_i,y_i\})dx^\prime

Parameters:
-----------
full_grid: numpy.ndarray
    array of points over which the smooth interpolation function is generated
nodes: numpy.ndarray
    array of nodes (x_i)
parameters: numpy.ndarray
    array of parameters corresponding with the nodes (y_i)
gamma: float
    smoothing parameter

Returns:
--------
smooth_interpolation: numpy.ndarray
    smooth_interpolation array that corresponds to the full_grid input
