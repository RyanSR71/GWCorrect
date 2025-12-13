Welcome to GWCorrect's documentation!
=====================================

Package version:
|version|

This package provides infrastructure and useful functions for modifying gravitational waveforms. With this package, you can correct waveforms for waveform uncertainty during parameter estimation runs (`Read 2023 <https://arxiv.org/abs/2301.06630v2>`_) and apply an inspiral phase correction to search for theories of gravity beyond general relativity (`Bonilla et. al. 2023 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.024015>`_). This work is in collaboration with and derived from past work by Dr. Jocelyn Read and Dr. Marceline Bonilla. 

This package relies on: `bilby <https://bilby-dev.github.io/bilby/index.html>`_, `PESummary <https://lscsoft.docs.ligo.org/pesummary/stable/index.html>`_, `pycbc <https://pycbc.org/pycbc/latest/html/>`_, `lal <https://lscsoft.docs.ligo.org/lalsuite/lal/modules.html>`_, `scipy <https://scipy.org/>`_, `numpy <https://numpy.org/>`_, `matplotlib <https://matplotlib.org/stable/index.html>`_, and `tqdm <https://tqdm.github.io/>`_.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation

.. toctree::
   :maxdepth: 1
   :caption: Tutorials:

   notebooks/Parameterization_and_Waveform_Differences.ipynb
   notebooks/ppE_Tutorial.ipynb

.. toctree::
   :maxdepth: 1
   :caption: API:

   wfu/index
   ppE/index

.. meta::
   <meta name="google-site-verification" content="IgwxvFJBW4t67HMJWbM1bx2dCEtdn8_lSaVh2kR8PIs" />
