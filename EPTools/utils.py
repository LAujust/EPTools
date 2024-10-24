import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as c
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from pyasassn.client import SkyPatrolClient
from astropy.io.votable import parse_single_table
from astropy.utils.data import conf

def keV2Hz(kevs):
    return 2.417990504024e+17 * kevs

def keV2T(kevs):
    return 11604525.00617 * kevs

def lam2Hz(lams):
    'lams[float/np.array]:  wavelength in Angstrom'
    lams = lams * u.Angstrom
    return (c.c/lams).cgs.value

def Hz2lam(nu):
    'Output:    wavelength in Angstrom'
    nu = nu * u.Hz
    return (c.c/nu).cgs.value * 1e8

def keV2lam(kevs):
    return Hz2lam(keV2Hz(kevs))

def mag2flx(mags,lam_ref=None,FWHM=None):
    """
    Output[flx]:    in erg s^-1 cm^-2 Hz^-1
    Output[flx]:    in erg s^-1 cm^-2 if wavelength and FWHM are given
    """
    if lam_ref is not None and FWHM is not None:
        delta_nu = lam2Hz(lam_ref-FWHM) - lam2Hz(lam_ref+FWHM)
        return delta_nu * 10**(-0.4*(mags+48.6))
    else:
        return 10**(-0.4*(mags+48.6))

def flx2mag(flxs):
    """
    flxs[float/array]:  in erg s^-1 cm^-2 Hz^-1
    """
    return -2.5*np.log10(flxs) - 48.6

def flx2lum(f,d):
    """
    f[flux]:        erg/s/cm^2
    d[distance]:    pc  
    output[L]:      erg/s
    """
    f = f * u.erg/u.s/(u.cm)**2
    d = d * u.pc
    L = 4*pi*d**2*f
    return L.cgs.value

def lum2flux(L,d):
    """
    L[luminosity]:  erg/s
    d[distance]:    pc
    output[f]:      erg/s/cm^2
    """
    L = L * u.erg/u.s
    f = f * u.erg/u.s/(u.cm)**2
    d = d * u.pc
    f = L/(4*pi*d**2)
    return f.cgs.value 

def data2acs(data,out_dir):
    t, cr, cr_err = data[:,0], data[:,1], data[:,2]
    
    tstart = t
    tstop = t + (t[1]-t[0])
    
    third = cr*(tstop-tstart)
    fourth = cr_err*(tstop-tstart)
    
    out = np.vstack((tstart,tstop,third,fourth)).T
    
    print(out)
    np.savetxt(prefix+'.txt',out)
    
    
    
#========================================================================#

pi = np.pi
band_wavelength = {
    'R':(6427,1298),
    'V':(5345,840),
    'g':(4782,1419),
    'r':(6217,1327),
    'i':(7532,1244),
    'z':(8685,1024),
    'L':(5523,2330),
}