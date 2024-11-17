import numpy as np
import pandas as pd
import glob
import sys, os
import astropy.units as u
import astropy.constants as c
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from pyasassn.client import SkyPatrolClient
from astropy.io.votable import parse_single_table
from astropy.utils.data import conf
import sncosmo
from ligo.gracedb.rest import GraceDb
import ligo.skymap
import simsurvey

sys.path.append('$HEADAS/lib/python')

def keV2Hz(kevs):
    return 2.417990504024e+17 * kevs

def Hz2keV(nu):
    return nu/2.417990504024e+17

def keV2T(kevs):
    return 11604525.00617 * kevs

def keV2erg(kevs):
    return 1.602e-9 * kevs


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

def lam2keV(lams):
    return(Hz2keV(lam2Hz(lams)))

def mag2flx(mags,lam_ref=None,FWHM=None):
    """
    Output[flx]:    in erg s^-1 cm^-2 Hz^-1
    Output[flx]:    in erg s^-1 cm^-2 if wavelength and FWHM are given
    """
    if lam_ref is not None and FWHM is not None:
        delta_nu = lam2Hz(lam_ref-0.5*FWHM) - lam2Hz(lam_ref+0.5*FWHM)
        return delta_nu * 10**(-0.4*(mags+48.6)), delta_nu
    else:
        return 10**(-0.4*(mags+48.6))
    
def mag2flx_sncosmo(mag,band):
    ab = sncosmo.get_magsystem('ab')
    flx = ab.band_mag_to_flux(mag,band)
    #if you want to convert to erg you should multiple by e=(c.h*nu*u.Hz).cgs
    return flx # in photons / s / cm^2


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
    #f = f * u.erg/u.s/(u.cm)**2
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
    np.savetxt(out_dir,out)


def retrive_gracedb(query=''):
    far = 1 #per yr
    far = far/31557600
    far = 'far < 3.17e-8'
    client = GraceDb()
    events = client.superevents(query=query)
    event_messages = {}
    for event in events:
        #print(f"Event ID: {event['superevent_id']}")
        id = event['superevent_id']
        try:
            circular = client.files(id,id+'-update.json').json()
        except:
            circular = 'No update circular'
        event_messages[id] = circular

    return event_messages

def X2mAB(F_X,E_min,E_max):
    """
    Convert observed X-ray flux to AB magnitude.
    
    Parameters:
    - F_X : float : Observed X-ray flux in erg/s/cm² (integrated over the band)
    - E_min : float : Minimum energy of the X-ray band in keV
    - E_max : float : Maximum energy of the X-ray band in keV
    
    Returns:
    - m_AB : float : AB magnitude corresponding to the observed X-ray flux
    """

    # Constants
    h = 6.626e-27  # Planck's constant in erg·s
    erg_to_keV = 1.60218e-9  # Conversion factor from erg to keV
    ab_zero_flux_density = 3.63e-20  # Zero-point for AB magnitude in erg/s/cm²/Hz

    # Step 1: Calculate the central energy (E_c) and bandpass (ΔE) in keV
    E_c = (E_min + E_max) / 2  # Central energy in keV
    delta_E = E_max - E_min    # Energy range in keV

    # Step 2: Convert central energy to frequency (ν_c)
    nu_c = (E_c * erg_to_keV) / h  # Frequency in Hz

    # Step 3: Convert integrated flux (F_X) to flux density per unit frequency (f_ν)
    f_nu = (F_X / delta_E) * (h / (E_c * erg_to_keV) ** 2)  # Flux density in erg/s/cm²/Hz

    # Step 4: Convert flux density to AB magnitude (m_AB)
    m_AB = -2.5 * np.log10(f_nu / ab_zero_flux_density)

    return m_AB


def EPexpo2mAB(expo,instrument:str):
    """
    expo[float, array]:     exposure time in second
    instrument[str]:        WXT or FXT
    """
    if instrument == 'WXT':
        F_X = 6.13670736e-09 * expo**-0.833761923
        E_min, E_max = 0.5, 4
    elif instrument == 'FXT':
        F_X = 5.49101782e-11 * expo**-0.812917739
        E_min, E_max = 0.5, 10
    else:
        raise KeyError('Please type correct instrument name!')
    
    return X2mAB(F_X,E_min,E_max)
    

def BNS_ejecta_mass(Mc, q, R1, R2, Mtov):
    pass

def NSBH_ejecta_mass(M_BH, M_NS, Chi, R_NS):
    pass

    
    
    
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
    'u':(3467,668)
}

x_instrument_energy_range = {
    'EP-WXT':[0.5,4],
    'EP-FXT':[0.4,10],
    'XRT':[0.4,10]
}