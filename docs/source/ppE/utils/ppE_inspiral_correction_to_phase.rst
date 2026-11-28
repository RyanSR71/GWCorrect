ppe_inspiral_correction_to_phase
================================
.. code-block:: python

  GWCorrect.ppE.utils.ppe_inspiral_correction_to_phase(freqs, total_mass, ppe_beta, ppe_b)

Computes the ppE phase correction for the inspiral regime Mf<Mf_IM. This function is defined assuming geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses. Likewise for if the total mass is given seconds.

.. math::

  \Delta\phi=\beta u^b

Parameters:
-----------
freqs: numpy.ndarray
  array of frequencies in geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses, likewise for if the total mass is given seconds
total_mass: float
  total binary mass in geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses, likewise for if the total mass is given seconds
ppe_beta: float
  ppE phase correction amplitude parameter
ppe_b: float
  ppE PN order parameter

Returns:
--------
phase_correction: numpy.ndarray
  ppE inspiral phase correction
