TFDG
====

.. code-block:: python

  class GWCorrect.wfu.prior.TFDG(mu_1,mu_2,sigma_1,sigma_2,minimum,maximum,name=None,latex_label=None)

Truncated Flattened Double Gaussian Prior: a double Gaussian prior that is uniform between the two means 
  - alternate prior for xi_0 if mu_down is greater than xi_min

.. math::

  \mathrm{TFDG}(x)=\begin{cases}
        N\exp\left(-\frac{1}{2}\left(\frac{x-\mu_1}{\sigma_1}\right)^2\right) & x\in[a,\mu_1] \\
        N & x\in(\mu_1,\mu_2) \\
        N\exp\left(-\frac{1}{2}\left(\frac{x-\mu_2}{\sigma_2}\right)^2\right) & x\in[\mu_2,b] 
    \end{cases}

.. math::

  N =\bigg(\sqrt{\frac{\pi}{2}}\bigg(\sigma_2\mathrm{erf}\left[\frac{b-\mu_2}{\sqrt{2}\sigma_2}\right]-\sigma_1\mathrm{erf}\left[\frac{a-\mu_1}{\sqrt{2}\sigma_1}\right]\bigg)+\mu_2-\mu_1\bigg)^{-1}

.. math::

  X_\mathrm{TFDG}(p)=\begin{cases}
        \sqrt{2}\sigma_1\mathrm{erf}^{-1}\left[\sqrt{\frac{2}{\pi}}\frac{p-A_1}{ N\sigma_1}\right]+\mu_1 & p\in[0,A_1] \\
        (p+N\mu_1-A_1)/N & p\in(A_1,A_1+A_2) \\
        \mu_2+\sqrt{2}\sigma_2\mathrm{erf}^{-1}\left[\sqrt{\frac{2}{\pi}}\frac{p-A_1-A_2}{ N\sigma_2}\right] & p\in[A_1+A_2,1]
    \end{cases}

.. math::

  A_1=\sqrt{\frac{\pi}{2}}N\sigma_1\mathrm{erf}\left[\frac{\mu_1-a}{\sqrt{2}\sigma_1}\right]

.. math::

  A_2=N(\mu_2-\mu_1)

Parameters:
-----------
mu_1: float
    mean for the first Gaussian
sigma_1: float
    standard deviation for the first Gaussian
mu_2: float
    mean for the second Gaussian
sigma_2: float
    standard deviation for the second Gaussian
minimum: float
    minimum value
maximum: float
    maximum value
name: string, optional, (None)
    name of the prior
latex_label: string, optional, (None)
    latex label of the prior

Returns:
--------
TFDG prior object that is compatible with bilby.
