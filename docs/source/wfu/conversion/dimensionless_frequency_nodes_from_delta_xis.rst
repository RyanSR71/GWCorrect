dimensionless_frequency_nodes_from_delta_xis
============================================

.. code-block:: python

  GWCorrect.wfu.conversion.dimensionless_frequency_nodes_from_delta_xis(xi_0,xi_n,delta_xis)

Computes the dimensionless frequency nodes from the outer frequency nodes and inner frequency node parameters.

.. math::

  \Xi_{k\ell}=\delta_{k\ell}-\delta_{k,l-1}\delta\xi_k-\delta_{k,l+1}(1-\delta\xi_k)

.. math::

  C_\ell=\delta_{\ell 1}(1-\delta\xi_\ell)\ln(\xi_0)+\delta_{\ell,n-1}\delta\xi_\ell\ln(\xi_n)

.. math::

  \xi_k=\exp\left([\Xi^{-1}]_{k\ell}C^\ell\right)

Parameters:
-----------
xi_0: float
    position of zeroth frequency node
xi_n: float
    position of nth frequecy node
delta_xis: list
    list of inner frequency node parameters

Returns:
--------
dimensionless_frequency_nodes: numpy.ndarray
    array of dimensionless frequency nodes from xi_0 to xi_n
