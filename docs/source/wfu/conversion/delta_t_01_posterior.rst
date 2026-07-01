delta_t_01_posterior
====================

.. code-block:: python

  GWCorrect.wfu.utils.delta_t_01_posterior(result,n,xi_max)

Computes delta_t_01 from a bilby result object. Used to add delta_t_01 to the posterior.

.. math::

  \Delta t_{01} = \frac{nGM}{c^3(\xi_\mathrm{max}-\xi_0)\delta\tilde\xi}

Parameters:
-----------
result: bilby.core.result.Result
    bilby result object
n: int
    number of frequency nodes excluding 0
xi_max: float
    upper bound on the dimensionless frequency band

Returns:
--------
delta_t_01_posterior: numpy.ndarray
    delta_t_01 posterior samples
