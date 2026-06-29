gaussian_parameters_from_A_ASD_solutions
========================================

.. code-block:: python

  GWCorrect.wfu.utils.gaussian_parameters_from_A_ASD_solutions(lower_xis,upper_xis,xi_max)

Computes the gaussian fits for the xi_0 and delta_xi_tilde priors. To be used in conjunction with GWCorrect.wfu.utils.A_ASD_solutions

Parameters:
-----------
lower_xis: list
    lower_xis from GWCorrect.wfu.utils.A_ASD_solutions
upper_xis: list
    upper_xis from GWCorrect.wfu.utils.A_ASD_solutions
xi_max: float
    upper bound on the dimensionless frequency band

Returns:
--------
mu1: float
    mean of the lower_xis distribution
sigma1: float
    standard deviation of the lower_xis distribution
mu2: float
    mean of the upper_xis distribution
sigma2: float
    standard deviation of the upper_xis distribution
mu: float
    mean of the delta_xi_tilde distribution
sigma: float
    standard deviation of the delta_xi_tilde distribution
