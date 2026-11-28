apply_ppe_correction
====================
.. code-block:: python

  GWCorrect.ppE.utils.apply_ppe_correction(strain, frequency_array, total_mass, ppe_beta, ppe_b, 
  ppe_delta_epsilon=0.0, Mfreq_IM = 0.018, ratio_f_MR_to_f_RD=0.75, aligned=True)

Applies the parameterized post-Einsteinian corrections to a waveform in the frequency domain such that the phase is first-differentiable across the different frequency regimes, where a different correction is applied in each one.

.. math::

  \Delta\phi_\mathrm{ppE}=\begin{cases} \Delta\phi_\mathrm{ppE}^\mathrm{IM} & f<f_\mathrm{IM}\\ \Delta\phi_\mathrm{Int} & f_\mathrm{IM}\leq f\leq f_\mathrm{MR}\\ \Delta\phi_\mathrm{ppE}^\epsilon & f > f_\mathrm{MR} \end{cases}

Parameters:
-----------
strain: numpy.ndarray
  Gravitational-wave strain
frequency_array: numpy.ndarray
  frequency array that goes with the strain
total_mass: float
  The total mass of the binary, in solar masses.
ppe_beta: float
  ppE phase correction amplitude
ppe_b: float
  ppE PN order parameter
ppe_delta_epsilon: float, optional, (0.0)
  ppE ringdown correction parameter
Mfreq_IM: float, optional, (0.018)
  The dimensionless frequency at which the inspiral regime ends.
ratio_f_MR_to_f_RD : float, optional, (0.75)
  The fraction of the ringdown frequency at which the intermediate regime ends and becomes the merger-ringdown regime.
aligned: bool, optional, (True)
  If True, the output waveform will be aligned to the input waveform in time and phase
    
    
Returns:
--------
new_fd_waveform : FrequencySeries
  The frequency domain waveform of the coalescence, after applying the ppE corrections.
