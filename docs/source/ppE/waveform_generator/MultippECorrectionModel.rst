MultippECorrectionModel
==================

.. code-block:: python

   GWCorrect.ppE.waveform_generator.MultippECorrectionModel(parameters,**kwargs)

EXPERIMENTAL: Modified lal binary black hole frequency domain source model to apply multiple ppE corrections to one waveform

.. math::

   \mu_{n\mathrm{ppE}}(f)=\mu(f;\vartheta)\exp\left(i\sum_{k=1}^n\Delta\phi_\mathrm{ppE}^\epsilon(f;\tilde\beta_k,\delta\tilde\epsilon_k,b_k)\right)

Make sure to include all "beta_tilde_k", "delta_epsilon_tilde_k", and all "b_k" in the prior.

New Parameters:
---------------
beta_tildes: list
  list of beta_tilde parameters
delta_epsilon_tildes: list   
  list of delta_epsilon_tilde parameters
bs: list
  list of b values
Mfreq_IM: float, optional, (0.018)
  The dimensionless frequency at which the inspiral regime ends.
ratio_f_MR_to_f_RD : float, optional, (0.75)
  The fraction of the ringdown frequency at which the intermediate regime ends and becomes the merger-ringdown regime.
aligned: bool, optional, (True)
  if True, ppE phase correction will not include time or phase shifts relative to the GR waveform


Returns:
--------
model_strain: dict
  plus and cross polarizations of the frequency domain strain
