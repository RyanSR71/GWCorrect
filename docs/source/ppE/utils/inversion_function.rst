inversion_function
==================
.. code-block:: python

  GWCorrect.ppE.utils.inversion_function(a,b,x)

Function that returns 1 when x is outside [a,b] and -1 when within [a,b]

.. math::

  \mathcal{I}_a^b(x)=\begin{cases} 1 & x<a\\ -1 & a\leq x\leq b \\ 1 & x>b \end{cases}

Parameters:
-----------
a: float
  lower bound on inversion interval
b: float
  upper bound on inversion interval
x: float
  input value

Returns:
--------
I: int
  1 or -1 depending on whether or not x was in [a,b]
