Quick Start
===========

A quick guide for using EPTools as a Transient Advocate (TA) and for subsequent data analysis.

Installation
------------

Prerequisites
~~~~~~~~~~~~~

- Python >= 3.8
- `Heasoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ with PyXspec correctly configured in your ``$PYTHONPATH``

From Source
~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/LAujust/EPTools.git
    cd EPTools
    pip install -r requirements.txt
    pip install .

The core dependencies (numpy, scipy, matplotlib, pandas, astropy, astroquery) are installed automatically.
For GW-related functionality, additional packages are listed in ``requirements.txt``.

WXT Lightcurve and Spectrum Analysis
-------------------------------------

Process WXT data and obtain lightcurve, spectrum, and spectral fitting for a target of opportunity (ToO):

.. code-block:: python

    import EPTools

    root = 'path/to/your/data'
    obsid = 'ep01709201259wxt37'
    snum = 's6'

    EPTools.TA_quick(
        obsid, snum, root,
        binsize=100,
        grp=True,
        group=1,
        ins='WXT',
        plotstyle='line',
        chatter=0,
    )

Lightcurve Extraction
---------------------

Use ``EPTools.analysis`` to extract source and background lightcurves from event files via ``xselect``:

.. code-block:: python

    from EPTools.analysis import extract_curve_sh, subtract_curve

    extract_curve_sh(
        evtfile='path/to/evt.fits',
        src_reg='src.reg',
        bkg_reg='bkg.reg',
        binsize=10,
    )

    subtract_curve('src.lc', 'bkg.lc', alpha=1/12.0, out_file='net.lc')

Spectral Fitting with Xspec
---------------------------

.. code-block:: python

    from EPTools.fit import xspec_fitting

    xspec_fitting(
        'src.pi', 'model_name',
        grp=True,
        rebin=2,
        stat='cstat',
        instrument='WXT',
        chatter=10,
    )

Plotting
--------

.. code-block:: python

    from EPTools.plot import lcurve_plot, xspec_plot

    # Plot background-subtracted lightcurve
    lcurve_plot('src.lc', 'bkg.lc', binsize=10, scale=1/12.0)

    # Plot Xspec fit results
    xspec_plot(data, save_dir='./plots', plotstyle='step')

Cross-matching
--------------

.. code-block:: python

    from EPTools.crossmatch import match_cat

    result = match_cat('source.cat', 'catalog.cat', radius=5.0)

Unit Conversions
----------------

.. code-block:: python

    from EPTools.utils import keV2erg, mag2flx, flx2lum

    # Convert keV to erg
    energy_erg = keV2erg([0.5, 2.0])

    # Convert AB magnitude to flux
    flux = mag2flx([20.0, 21.5])

    # Convert flux to luminosity
    lum = flx2lum(flux=1e-14, d=100)  # distance in Mpc
