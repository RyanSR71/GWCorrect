TotalMassConstraint
===================

.. code-block:: python

  GWCorrect.wfu.prior.TotalMassConstraint(*,name,minimum_frequency,maximum_frequency,
  latex_label=r'$M$',boundary=None,unit=r'$\mathrm{M}_\odot$',xi_min=0.018,xi_max=0.318)

Generates a bilby prior that constrains the total mass to ensure that the waveform correction is always within the frequency band.

.. math::

  \frac{c^3\xi_\mathrm{high}}{Gf_\mathrm{high}}<M<\frac{c^3\xi_\mathrm{low}}{Gf_\mathrm{low}}

Parameters:
-----------
name: string
  name of prior
minimum_frequency: float
  lower bound on frequency band in Hz
maximum_frequency: float
  upper bound on the frequency band in Hz
latex_label: string, optional, (r'$M$')
  label for the parameter in LaTeX
boundary: string, optional, (None)
  boundary condition for the prior
unit: string, optional, (r'$\mathrm{M}_\odot$')
  label for the unit of the parameter; default is solar mass
xi_min: float, optional, (0.018)
  lower bound on the dimensionless frequency band
xi_max: float, optional, (1/pi, 0.318...)
  upper bound on the dimensionless frequency band

Returns:
--------
total_mass_prior: bilby.core.prior.base.Constraint
  bilby constraint prior object for the total mass
