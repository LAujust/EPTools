import EPTools
from astropy.coordinates import SkyCoord
from EPTools.utils import *
import numpy as np
from astropy.cosmology import Planck18
import matplotlib.pyplot as plt
import sys, os
sys.path.append('$HEADAS/lib/python')
from xspec import *


"""
!!!!!
A problem is that the.arf and .rmf should be in the same dir of .py if you want 
to use .pi spectrum generated by ftgrppha
!!!!!
"""

root_dir = '/Users/liangrunduo/Desktop/Aujust/NAOC/EP/EPTools/data/'
s1 = Spectrum(root_dir+'PC.pi')
s1.ignore("0.0-0.5 10.0-**")
AllData.ignore("bad")

m = Model('tbabs*cflux*powerlaw')
m.setPars({1:0.05,6:1})
m.TBabs.nH.frozen = True
m.powerlaw.PhoIndex.frozen = True
Fit.perform()
Plot.device = '/null'
Plot('data')
# Get coordinates from plot:
chans = Plot.x()
rates = Plot.y()
folded = Plot.model()

# Plot using Matplotlib:
plt.plot(chans, rates, 'ro', chans, folded)
plt.xlabel('channels')
plt.ylabel('counts/cm^2/sec/chan')
plt.savefig('myplot')

