WaveformCorrectionModelBBH
==========================

.. code-block:: python

   GWCorrect.wfu.waveform_generator.WaveformCorrectionModelBBH(parameters,**kwargs)

Modified lal binary black hole frequency domain source model to include waveform error corrections in the strain calculation.

.. math::

   \mu_\mathcal{C}(f;\vartheta,\mathrm{A},\Phi)=\mu(f;\vartheta)(1+\Delta\mathcal{A}_\mathrm{SI}(f;\mathrm{A}))\exp(i\Delta\phi_\mathrm{SI}(f;\Phi))

Make sure to include "xi_0", "xi_n", and all necessary "alpha_tilde", "phi_tilde', and "delta_xi" parameters in the prior.

New Parameters:
---------------
xi_0: float
   position of the zeroth frequency node in dimensionless frequency
xi_n: float
   position of the nth frequency node in dimensionless frequency
delta_xis: list
   list of inner frequency node parameters
alpha_tildes: list
   list of rescaled amplitude correction parameters
phi_tildes: list
   list of rescaled phase correction parameters
sigma_A_spline: scipy.interpolate.CubicSpline
   spline object from which to calculate the sigma values for the amplitude correction; put this in the waveform arguments dictionary
sigma_phi_spline: scipy.interpolate.CubicSpline
   spline object from which to calculate the sigma values for the phase correction; put this in the waveform arguments dictionary
xi_max: float, optional, (1/pi)
   absolute upper bound on dimensionless frequency

Returns:
--------
model_strain: dict
   plus and cross polarizations of the frequency domain strain
