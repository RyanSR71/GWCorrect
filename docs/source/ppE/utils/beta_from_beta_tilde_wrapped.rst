beta_from_beta_tilde_wrapped
============================
.. code-block:: python

  GWCorrect.ppE.utils.beta_from_beta_tilde_wrapped(beta_tilde, f_min_Hz, Mf_max, exponent_b, Mf_IM, total_mass)

Converts between beta_tilde and beta.

.. math::

  \beta=\tilde{\beta}u_\mathrm{IM}^{-b}\left(\Phi_b-\frac{1}{2}\vec{p_b}^\mathbf{T} P^{-1}\vec{p_b}\right)^{-\frac{1}{2}}


Parameters:
-----------
beta_tilde : float
  The rescaled beta term (which can be computed from a desired mismatch).
f_min_Hz : float
  The lower frequency cutoff of the waveform, in units of Hz.
Mf_max : float
  The upper frequency cutoff of the waveform, in dimensionless units.
exponent_b : float
  The ppE exponent b.
Mf_IM: float
  The dimensionless frequency at which the inspiral regime ends.
total_mass: float
  The total mass of the binary in units of solar masses.

                           
Returns:
--------
beta : float
  The ppE parameter beta.
