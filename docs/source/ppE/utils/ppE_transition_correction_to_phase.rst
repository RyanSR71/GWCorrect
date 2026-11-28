ppe_transition_correction_to_phase
==================================
.. code-block:: python

  GWCorrect.ppE.utils.ppe_transition_correction_to_phase(freqs, phi_IM, phi_MR, dphi_IM, dphi_MR, f_IM, f_MR)

Computes the post-inspiral ppE phase correction. This function is defined assuming geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses. Likewise for if the total mass is given seconds.

.. math::

  \Delta\phi_\mathrm{Int}(f;\beta,\delta\epsilon)

Parameters:
-----------
freqs: numpy.ndarray
  array of frequencies in geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses, likewise for if the total mass is given seconds
phi_IM: float
  phase at Mf_IM
phi_MR: float
  phase at Mf_MR
dphi_IM: float
  ppE phase correction at Mf_IM
dphi_MR: float
  ppE phase correction at Mf_MR
f_IM: float
  the dimensionless frequency at which the inspiral regime ends
f_MR: float
  the dimensionless frequency at which the merger-ringdown regime starts

Returns:
--------
phase_correction: numpy.ndarray
  transition/intermediate ppE phase correction
