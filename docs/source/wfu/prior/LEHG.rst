LEHG
====

.. code-block:: python

  class GWCorrect.wfu.prior.LEHG(mu,sigma,minimum,maximum,name=None,latex_label=None)

Left Extended Half Gaussian Prior: a truncated Gaussian prior that is uniform between the minimum and the mean
  - alternate prior for xi_0 if mu_down is less than xi_min
  - alternate prior for delta_xi_tilde

.. math::

  \mathrm{LEHG}(x)=\begin{cases}
        N & x\in[a,\mu] \\
        N\exp\left(-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2\right) & x\in[\mu,b]
    \end{cases}

.. math::

  N=\left(\mu-a-\sqrt{\frac{\pi}{2}}\sigma\mathrm{erf}\left[\frac{\mu-b}{\sqrt{2}\sigma}\right]\right)^{-1}

.. math::

  X_\mathrm{LEHG}(p)=\begin{cases}
        p/N & p\in[0,N\mu) \\
        \mu+\sqrt{2}\sigma\mathrm{erf}^{-1}\left[\sqrt{\frac{2}{\pi}}\frac{p-N\mu}{ N\sigma}\right] & p\in[N\mu,1]
    \end{cases}

Parameters:
-----------
mu: float
    mean
sigma: float
    standard deviation
minimum: float
    minimum value
maximum: float
    maximum value
name: string, optional, (None)
    name of the prior
latex_label: string, optional, (None)
    latex_label of the prior

Returns:
--------
LEHG prior object that is compatible with bilby.
