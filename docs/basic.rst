Basic usage
====================

Calculate Luminosity
-----------------------------------------------

.. code-block:: python

    f = 1e-11 #erg/cm^2/s
    d = 500 #pc
    L = EPTools.flx2lum(f,d)
    print(f"L = {L:.2e} erg/s")


Estimate Significance and Upper Limits
------------------------------------------------

.. code-block:: python

    N_on = 10
    N_off = 30
    alpha = 1./12 #The ratio of size of source and background
    S = EPTools.li_ma_sigma(N_on, N_off, alpha)

.. code-block:: python

    Nsrc = 10 #photons within src region
    Nbkg = 30 #photons within bkg region
    factor = 2e-9 #ctr to flux conversion factor

    X_UL(Nsrc,Nbkg,exposure,alpha=1./12,factor=factor,CL = 0.9)


