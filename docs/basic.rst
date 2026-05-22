Basic Usage
===========

Calculate Luminosity
--------------------

.. code-block:: python

    import EPTools

    f = 1e-11  # erg/cm^2/s
    d = 500    # pc

    L = EPTools.flx2lum(f, d)
    print(f"L = {L:.2e} erg/s")

Estimate Significance and Upper Limits
--------------------------------------

Li & Ma significance:

.. code-block:: python

    import EPTools

    N_on = 10
    N_off = 30
    alpha = 1./12  # ratio of source to background region size
    S = EPTools.li_ma_sigma(N_on, N_off, alpha)
    print(f"Significance: {S:.2f} sigma")

Upper limit:

.. code-block:: python

    Nsrc = 10        # photons within source region
    Nbkg = 30        # photons within background region
    exposure = 1000  # seconds
    factor = 2e-9    # count-rate to flux conversion factor

    ul = EPTools.X_UL(Nsrc, Nbkg, exposure, alpha=1./12, factor=factor, CL=0.9)
    print(f"Upper limit: {ul:.2e} erg/cm^2/s")

BAT Instrument Response and Sensitivity
---------------------------------------

Calculate count-rate to flux conversion factor for a given spectrum:

.. code-block:: python

    from EPTools.fit import get_ctrt_to_flux

    factor = get_ctrt_to_flux(
        source_spec='powerlaw',
        energy_l=0.5,
        energy_h=4.0,
        nH_Galactic=0.05,
        PhoIndex=1.8,
        get_unabs=True,
        ins='WXT',
    )

Read and Plot Lightcurve
------------------------

.. code-block:: python

    from EPTools.utils import read_curve
    from EPTools.plot import lcurve_plot

    times, rates, errors = read_curve('src.lc', 'bkg.lc', binsize=10, scale=1./12)
    lcurve_plot('src.lc', 'bkg.lc', binsize=10, scale=1./12)
