ppECorrectionModel
==================

.. code-block:: python

   GWCorrect.ppE.waveform_generator.ppECorrectionModel(parameters,**kwargs)

Modified lal binary black hole frequency domain source model to include ppE correction in the strain calculation.

.. math::

   \mu_\mathrm{ppE}(f;\vartheta,\Phi)=\mu(f;\vartheta)\exp(i\Delta\phi_\mathrm{ppE}(f;\tilde\beta,\delta\tilde\epsilon,b))

Make sure to include "beta_tilde", "delta_epsilon_tilde", and "b" in the prior.

New Parameters:
---------------
beta_tilde: float
   rescaled ppE inspiral correction parameter
delta_epsilon_tilde: float
   rescaled ppE post-inspiral correction parameter
b: float
   PN order parameter
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
