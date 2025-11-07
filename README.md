# EPTools

[![Hits-of-Code](https://hitsofcode.com/github/LAujust/EPTOols?branch=master)](https://hitsofcode.com/github/LAujust/EPTOols/view?branch=master) 
[![License](https://img.shields.io/github/license/LAujust/EPTools.svg)](https://github.com/LAujust/EPTools/LICENSE)

An easy python package/code for Einstein Probe data analysis.

## Useful Links

[ASAS-SN](https://asas-sn.osu.edu/photometry) forced photometry.

[ATLAS](https://fallingstar-data.com/forcedphot/queue/) forced photometry. 

[Alerce](https://alerce.online/) and [Lasair](https://lasair-ztf.lsst.ac.uk/) are commonly used platform for ZTF sources. 

Heasarc tool provide a convenient [portal](https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/high_energy_source/high_energy_source.pl) to search history X-ray observations. You may also need check eROSITA via https://erosita.mpe.mpg.de/dr1/erodat/catalogue/search/.

A [tool](https://www.swift.ac.uk/analysis/nhtot/) used for obtain Galactic hydrogen column density $N_{\mathrm{H}}$.

[RapidGBM](https://hetools.xyz/RapidGBM) for check Fermi-GBM sky coverage. 

[DataFusion](https://nadc.china-vo.org/gwops/#/), a Multi-messenger tool for viewing EP instant coverage and ToO Plan. 

[BAT Target Search Portal](https://guano.swift.psu.edu/) offers a quick access to Swift-BAT instant coverage.

[Astro-Colibri](https://astro-colibri.com/) is a platform for viewing transient information of different type. 

[Legacy Survey](https://www.legacysurvey.org/): useful Gallery for checking deep optical image.

## Requirements

Heasoft should be installed, and the path of PyXspec need to be correctly added. 

## Quick Start

After installed `EPTools`, you can easily process WXT data by using `EPTools.TA_quick` function to obtain basic lightcurve, spectrum and fitting of specified source for TA:

```
import EPTools
root = 'path/to/your/data'
obsid = 'ep01709201259wxt37'
snum = 's6'
binsize = 100
grp = True
group = 1
nH = None
pha_file = 'PC.pi'  #if you have already grouped data
plotstyle = 'line'
ins = 'WXT'
chatter = 0

EPTools.utils.TA_quick(obsid,snum,root,grp=grp,pha_file=pha_file,group=group,nH=nH,
    binsize=binsize,rebin=rebin,rx=rx,plotstyle=plotstyle,chatter=chatter,ins=ins)

```

More informations can be found in the [documentation](https://eptools.readthedocs.io/en/latest/).

@copyright to Runduo Liang (Aujust)

