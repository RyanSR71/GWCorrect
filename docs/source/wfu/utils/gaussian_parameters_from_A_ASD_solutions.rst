gaussian_parameters_from_A_ASD_solutions
========================================

.. code-block:: python

  GWCorrect.wfu.utils.gaussian_parameters_from_A_ASD_solutions(lower_xis,upper_xis)

Computes the gaussian fits for the xi_0 and xi_n priors. To be used alongside GWCorrect.wfu.utils.A_ASD_solutions

Parameters:
-----------
lower_xis: list
    lower_xis from GWCorrect.wfu.utils.A_ASD_solutions
upper_xis: list
    upper_xis from GWCorrect.wfu.utils.A_ASD_solutions

Returns:
--------
mu_0: float
    mean of the lower_xis distribution
sigma_0: float
    standard deviation of the lower_xis distribution
mu_n: float
    mean of the upper_xis distribution
sigma_n: float
    standard deviation of the upper_xis distribution
