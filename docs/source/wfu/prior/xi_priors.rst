xi_priors
=========

.. code-block:: python

   GWCorrect.wfu.prior.xi_priors(waveform_generator,prior,psd_data,n,minimum_frequency,xi_min=0.018,
   xi_max=1/np.pi,xi_0_latex_label=r'$\xi_0$',delta_xi_tilde_latex_label=r'$\delta\tilde\xi$')

Generates xi_0 and delta_xi_tilde priors from a BBH prior and adds them to the original prior.

.. math::

   2\mathcal{A}(\xi;\vartheta)\sqrt{\frac{c^3\xi}{GM}}-\sqrt{S_n(\xi)}=0

.. math::

  \Pi(\xi_0)=\mathrm{TFDG}(\mu_1,\mu_2,\sigma_1,\sigma_2,\xi_\mathrm{min},\xi_\mathrm{max}),\ \Pi(\xi_0)=\mathrm{EHG}(\mu,\sigma,\xi_\mathrm{min},\xi_\mathrm{max})

.. math::

   \Pi(\delta\tilde\xi)=\mathrm{EHG}(\mu,\sigma,0,1)                                                                                             

Parameters:
-----------
waveform_generator: bilby.gw.WaveformGenerator
    bilby waveform generator object
prior: bilby.core.prior.PriorDict
    bilby prior dictionary
psd_data: numpy.ndarray
    array of power spectral density data; first column needs to be the frequency points and the second column needs to be the data
n: int
    number of frequency nodes
minimum_frequency: float
    lower bound on the frequency band (Hz)
xi_0_latex_label: string, optional, (r'$xi_0$')
    latex label for xi_0
delta_xi_tilde_latex_label: string, optional, (r'$delta tilde xi$')
    latex_label for delta_xi_tilde
xi_min: float, optional, (0.018)
    lower bound on the dimensionless frequency band
xi_max: float, optional, (1/pi)
    upper bound on the dimensionless frequency band
samples: int, optional, (1000)
    number of draws of amplitude to take to generate the priors

Returns:
--------
prior: bilby.core.prior.PriorDict
  input prior dictionary, but with the new xi_0 and delta_xi_tilde priors added
