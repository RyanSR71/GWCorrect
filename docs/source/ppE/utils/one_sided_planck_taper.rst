one_sided_planck_taper
======================
.. code-block:: python

  GWCorrect.ppE.utils.one_sided_planck_taper(e,n,N)

Computes a one-sided Planck taper window.

.. math::

  w(n)=\begin{cases} \left(1+\exp\left(\frac{1}{n}-\frac{1}{n-eN}\right)\right)^{-1} & n<eN\\ 1 & n\geq eN \end{cases}

Parameters:
-----------
e: float
  taper parameter
n: float
  input value
N: int
  number of samples
                           
Returns:
--------
w: float
  value of the Planck taper window at n
