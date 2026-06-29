REHG
====

.. code-block:: python

  class GWCorrect.wfu.prior.REHG(mu,sigma,minimum,maximum,name=None,latex_label=None)

Right Extended Half Gaussian Prior: a truncated Gaussian prior that is uniform between the mean and the maximum
  - alternate prior for delta_xi_tilde

.. math::

  \mathrm{REHG}(x)=\begin{cases}
        Ne^{-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2} & x\in[a,\mu]\\
        N & x\in(\mu,b]
    \end{cases}

.. math::

  N=\left(\sqrt{\frac{\pi}{2}}\sigma\mathrm{erf}\left[\frac{\mu-a}{\sqrt{2}\sigma}\right]+b-\mu\right)^{-1}

.. math::

  X_\mathrm{REHG}(p)=\begin{cases}
        \mu-\sqrt{2}\sigma\mathrm{erf}^{-1}\left[\mathrm{erf}\left(\frac{\mu-a}{\sqrt{2}\sigma}\right)-\sqrt{\frac{2}{\pi}}\frac{p}{N\sigma}\right] & p\in[0,A]\\
        \mu+(p-A)/N & p\in(A,1]
    \end{cases}

.. math::

  A=N\sigma\sqrt{\frac{\pi}{2}}\mathrm{erf}\left(\frac{\mu-a}{\sqrt{2}\sigma}\right)

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
REHG prior object that is compatible with bilby.
