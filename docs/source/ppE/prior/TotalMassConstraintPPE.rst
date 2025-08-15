TotalMassConstraintPPE
======================

.. code-block:: python

  GWCorrect.ppE.prior.TotalMassConstraintPPE(*,name,minimum_frequency,maximum_frequency,latex_label=r'$M$',
  boundary=None,unit=r'$\mathrm{M}_\odot$',Mf_IM=0.018)

Generates a bilby prior that constrains the total mass to ensure that the ppE correction is always within the frequency band.

.. math::

  M<\frac{c^3[Mf_\mathrm{IM}]}{Gf_\mathrm{min}}

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
Mf_IM: float, optional, (0.018)
  end of IMRPhenom inspiral regime in dimensionless units

Returns:
--------
total_mass_prior: bilby.core.prior.base.Constraint
  bilby constraint prior object for the total mass
