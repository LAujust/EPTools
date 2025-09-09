# EPTools

[![Hits-of-Code](https://hitsofcode.com/github/LAujust/EPTOols?branch=master)](https://hitsofcode.com/github/LAujust/EPTOols/view?branch=master) 
[![License](https://img.shields.io/github/license/LAujust/EPTools.svg)](https://github.com/LAujust/EPTools/LICENSE)

An easy python package/code for Einstein Probe data analysis.

## Requirements

Heasoft should be installed, and the path of PyXspec need to be correctly e

## Quick Start

After installed `EPTools`, you can easily process WXT data by using `EPTools.TA_quick` function to obtain basic lightcurve, spectrum and fitting of specified source:

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

@copyright to Runduo Liang (Aujust)

