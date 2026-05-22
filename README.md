# EPTools

[![Hits-of-Code](https://hitsofcode.com/github/LAujust/EPTOols?branch=master)](https://hitsofcode.com/github/LAujust/EPTOols/view?branch=master)
[![License](https://img.shields.io/github/license/LAujust/EPTools.svg)](https://github.com/LAujust/EPTools/LICENSE)

A Python package for Einstein Probe (EP) data analysis, also supporting other X-ray instruments (e.g., Swift-XRT).

## Installation

### Prerequisites

- Python >= 3.8
- [Heasoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/) with PyXspec correctly configured in your `$PYTHONPATH`

### Install via pip

```bash
git clone https://github.com/LAujust/EPTools.git
cd EPTools
pip install -r requirements.txt
pip install .
```

### Required Python packages

The core dependencies (numpy, scipy, matplotlib, pandas, astropy, astroquery) will be installed automatically. For GW-related functionality, additional packages are listed in `requirements.txt`.

## Quick Start

### WXT lightcurve and spectrum analysis

Process WXT data and obtain lightcurve, spectrum, and spectral fitting for a target of opportunity (ToO):

```python
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
```

### Lightcurve extraction

Use `EPTools.analysis` to extract source and background lightcurves from event files via `xselect`:

```python
from EPTools.analysis import extract_curve_sh, subtract_curve

extract_curve_sh(
    evtfile='path/to/evt.fits',
    src_reg='src.reg',
    bkg_reg='bkg.reg',
    binsize=10,
)

subtract_curve('src.lc', 'bkg.lc', alpha=1/12.0, out_file='net.lc')
```

### Spectral fitting with Xspec

```python
from EPTools.fit import xspec_fitting

xspec_fitting(
    'src.pi', 'model_name',
    grp=True,
    rebin=2,
    stat='cstat',
    instrument='WXT',
    chatter=10,
)
```

### Plotting

```python
from EPTools.plot import lcurve_plot, xspec_plot

# Plot background-subtracted lightcurve
lcurve_plot('src.lc', 'bkg.lc', binsize=10, scale=1/12.0)

# Plot Xspec fit results
xspec_plot(data, save_dir='./plots', plotstyle='step')
```

### Cross-matching

```python
from EPTools.crossmatch import match_cat

result = match_cat('source.cat', 'catalog.cat', radius=5.0)
```

### Unit conversions

```python
from EPTools.utils import keV2erg, mag2flx, flx2lum

# Convert keV to erg
energy_erg = keV2erg([0.5, 2.0])

# Convert AB magnitude to flux
flux = mag2flx([20.0, 21.5])

# Convert flux to luminosity
lum = flx2lum(flux=1e-14, d=100)  # distance in Mpc
```

## Modules

| Module | Description |
|--------|-------------|
| `EPTools.utils` | General utilities, unit conversions, TA quick-look, GCN retrieval |
| `EPTools.analysis` | Lightcurve/spectrum extraction via xselect, response file generation |
| `EPTools.fit` | Xspec fitting, dynesty Bayesian fitting |
| `EPTools.plot` | Lightcurve, spectrum, and FITS image plotting |
| `EPTools.crossmatch` | Source catalog cross-matching |
| `EPTools.gw_em` | Gravitational-wave EM counterpart tools, EP observation planning |
| `EPTools.model` | Physical models |

## Useful Links

- [ASAS-SN](https://asas-sn.osu.edu/photometry) forced photometry
- [ATLAS](https://fallingstar-data.com/forcedphot/queue/) forced photometry
- [Alerce](https://alerce.online/) and [Lasair](https://lasair-ztf.lsst.ac.uk/) for ZTF sources
- [HEASARC High-Energy Source Portal](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/high_energy_source/high_energy_source.pl) — historical X-ray observations
- [eROSITA DR1 Catalogue](https://erosita.mpe.mpg.de/dr1/erodat/catalogue/search/)
- [NH Column Density Tool](https://www.swift.ac.uk/analysis/nhtot/)
- [RapidGBM](https://hetools.xyz/RapidGBM) — Fermi-GBM sky coverage
- [DataFusion](https://nadc.china-vo.org/gwops/#/) — EP coverage and ToO planning
- [BAT Target Search Portal](https://guano.swift.psu.edu/) — Swift-BAT coverage
- [Astro-Colibri](https://astro-colibri.com/) — transient information aggregator
- [Legacy Survey](https://www.legacysurvey.org/) — deep optical imaging

## Documentation

Full documentation is available at [eptools.readthedocs.io](https://eptools.readthedocs.io/en/latest/).

---

&copy; Runduo Liang (Aujust)
