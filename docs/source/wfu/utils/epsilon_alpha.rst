epsilon_alpha
=============

.. code-block:: python

  GWCorrect.wfu.utils.epsilon_alpha(dAs,xi_0,delta_xi_tilde,xi_max,sigma_dA_spline)

Computes epsilon_alpha up to 1.

.. math::

  \epsilon_\alpha=n-\sum_{k=1}^n\theta(\alpha_k+1)

Parameters:
-----------
dAs: list or numpy.ndarray
    list of dA parameters including dA_0
xi_0: float
    starting frequency node parameter
delta_xi_tilde: float
    frequency node spacing parameter
xi_max: float
    upper bound on the dimensionless frequency band
sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
    scipy cubic spline object encoding the standard deviation of the dA prior

Returns:
--------
epsilon_alpha: float
    amplitude error parameter
