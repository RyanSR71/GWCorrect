dimensionless_frequency_nodes_from_xi_0_delta_xi_tilde
======================================================

.. code-block:: python

  GWCorrect.wfu.conversion.dimensionless_frequency_nodes_from_xi_0_delta_xi_tilde(xi_0,delta_xi_tilde,n,xi_max)

Converts xi_0 and delta_xi_tilde to frequency nodes in mass-weighted dimensionless frequency.

.. math::

  \xi_k=\xi_0\left(1+\left[\frac{\xi_\mathrm{max}-\xi_0}{\xi_0}\right]\delta\tilde\xi\right)^{k/n}

Parameters:
-----------
xi_0: float
    starting frequency node parameter
delta_xi_tilde: float
    frequency_node_spacing_parameter
n: int
    number of frequency nodes excluding 0
xi_max: float
    upper bound on the dimensionless frequency band

Returns:
--------
dimensionless_frequency_nodes: numpy.ndarray
    array of dimensionless frequency nodes from xi_0 to xi_n
