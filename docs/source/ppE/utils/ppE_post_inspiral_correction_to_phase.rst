ppe_post_inspiral_correction_to_phase
=====================================
.. code-block:: python

  GWCorrect.ppE.utils.ppe_post_inspiral_correction_to_phase(freqs, total_mass, ppe_beta, ppe_b, 
  ppe_delta_epsilon=0.0, Mfreq_IM=0.018)

Computes the post-inspiral ppE phase correction. This function is defined assuming geometrized units; if the total mass is in solar masses, the frequency must be in inverse solar masses. Likewise for if the total mass is given seconds.

.. math::

  \Delta\phi_\mathrm{ppE}^\epsilon(f\geq f_\mathrm{MR})=\beta u_\mathrm{IM}^b\left(1+\frac{b(\beta+\delta\epsilon)}{3\beta}(p^3-1)\right)

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
ppe_delta_epsilon: float, optional, (0.0)
  ppE merger correction parameter
Mf_IM: float, optional, (0.018)
  the dimensionless frequency at which the inspiral regime ends

Returns:
--------
phase_correction: numpy.ndarray
  post-inspiral ppE phase correction
