conversion
==========

.. code-block:: python

  GWCorrect.wfu.prior.conversion(parameters,xi_max,n)

Computes the total mass and delta_t_01 from a set of system parameters. Used with the bilby.core.prior.PriorDict() function

.. math::

  \Delta t_{01}\equiv\frac{1}{\mathrm{f}_1-\mathrm{f}_0}

Parameters:
-----------
parameters: dictionary
  dictionary with binary black hole parameters
xi_max: float
  absolute upper bound on dimensionless frequency
n: int
  number of waveform correction parameters

Returns:
--------
parameters: dictionary
  input dictionary, but with the total system mass and delta_t_01 added
