import numpy as np
import astropy.units as u
import astropy.constants as c
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from pyasassn.client import SkyPatrolClient
from astropy.io.votable import parse_single_table

pi = np.pi
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