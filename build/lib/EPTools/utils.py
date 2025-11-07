import numpy as np
import pandas as pd
import sys, os, glob, re
import astropy.units as u
import subprocess
import astropy.constants as c
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table
from astropy.utils.data import conf
from ligo.gracedb.rest import GraceDb
import ligo.skymap
from scipy.special import gammainc, gammaincc, gammaincinv
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
sys.path.append('/Users/liangrunduo/heasoft-6.34/aarch64-apple-darwin23.5.0/lib/python')
sys.path.append('$HEADAS/lib/python')
import xspec as xs # type: ignore
from .plot import *
from .fit import *



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
    
# def mag2flx_sncosmo(mag,band):
#     '''
#     mag[float]:  AB magnitude
#     band[str]:   band name customized to sncosmo
    
#     Output:
#     flux[float]:   in erg/s/cm^2
#     '''
#     ab = sncosmo.get_magsystem('ab')
#     flx = ab.band_mag_to_flux(mag,band)
#     band = sncosmo.get_bandpass(band)
#     wave_eff = band.wave_eff
#     nu_eff = lam2Hz(wave_eff)
#     e = (c.h*nu_eff*u.Hz).cgs.value
#     #if you want to convert to erg you should multiple by e=(c.h*nu*u.Hz).cgs
#     return e*flx


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

def lcurve2pha(data_dir,out_dir):
    data = np.loadtxt(data_dir,skiprows=3)
    t, t_err, cr, cr_err = data[:,0], data[:,1], data[:,2], data[:,3]
    
    tstart = t - t_err[0]
    tstop = t + t_err[0]
    
    third = cr*(tstop-tstart)
    fourth = cr_err*(tstop-tstart)
    
    out = np.vstack((tstart,tstop,third,fourth)).T
    np.savetxt(out_dir,out)
    return out

def fplot2pha(data_dir,out_dir):
    data = np.loadtxt(data_dir,skiprows=3)
    t, cr, cr_err = data[:,0], data[:,1], data[:,2]
    
    tstart = t
    tstop = t + (t[1]-t[0])
    
    third = cr*(tstop-tstart)
    fourth = cr_err*(tstop-tstart)
    
    out = np.vstack((tstart,tstop,third,fourth)).T
    
    np.savetxt(out_dir,out)
    return out

def li_ma_sigma(N_on, N_off, alpha):
    if N_on <= 0 or N_off <= 0:
        return 0.0  # or np.nan
    term1 = N_on * np.log((1 + alpha) / alpha * (N_on / (N_on + N_off)))
    term2 = N_off * np.log((1 + alpha) * (N_off / (N_on + N_off)))
    return np.sqrt(2 * (term1 + term2))

def X_UL(Nsrc,Nbkg,exposure,alpha=1/12,factor=1e-9,CL = 0.9):
    B = Nbkg * alpha
    C_1 = gammaincc(Nsrc+1,B)
    sC_1 = gammainc(Nsrc+1,B)
    part1 = gammaincinv(Nsrc+1,CL*C_1+sC_1)
    UL = part1 - B
    return UL * factor / exposure
    
def get_ctrt_to_flux(source_spec, energy_l, energy_h, nH_Galactic, PhoIndex, get_unabs=True,nH_Intrinsic=0, ins='WXT'):
    """_summary_
    Args:
        source_spec (_type_): _description_
        energy_l (_type_): _description_
        energy_h (_type_): _description_
        nH (_type_): _description_
        PhoIndex (_type_): _description_
        get_unabs (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    xs.Xset.chatter = 0
    xs.Xset.logChatter = 0
    xs.Fit.statMethod = 'cstat'
    spec = xs.Spectrum(source_spec)
    spec.ignore('**-%.1f %.1f-**'%(energy_l, energy_h))
    spec.background = None
    model = xs.Model('tbabs*ztbabs*cflux*powerlaw')
    if not get_unabs:
        model = xs.Model('cflux*tbabs*ztbabs*powerlaw')
    model.TBabs.nH.values = nH_Galactic
    model.zTBabs.nH.values = nH_Intrinsic
    model.zTBabs.Redshift.values = 0
    model.cflux.Emin = energy_l
    model.cflux.Emax = energy_h
    model.cflux.lg10Flux = -9.0
    model.powerlaw.PhoIndex = PhoIndex
    model.powerlaw.norm = 1.0
    model.powerlaw.norm.frozen = True
    fakeit_kwargs = {}
    fakeit_kwargs['exposure'] = 10000
    fakeit_kwargs['correction'] = 1.0
    fakeit_kwargs['backExposure'] = 1.0
    fakeit_kwargs['fileName'] = 'temp_fake.pha'
    xs.AllData.fakeit(1, xs.FakeitSettings(**fakeit_kwargs))
    spec = xs.AllData(1)
    spec.ignore('**-%.1f %.1f-**'%(energy_l, energy_h))
    spec.show()
    model.show()
    ctrt = spec.rate[-1]
    flux = 10**(model.cflux.lg10Flux.values[0])
    xs.AllModels.clear()
    xs.AllData.clear()
    return flux/ctrt


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

def X2mAB(F_X,nu_min,nu_max):
    """
    Convert observed X-ray flux to AB magnitude.
    
    Parameters:
    - F_X : float : Observed X-ray flux in erg/s/cm² (integrated over the band)
    - E_min : float : Minimum energy of the X-ray band in keV
    - E_max : float : Maximum energy of the X-ray band in keV
    
    Returns:
    - m_AB : float : AB magnitude corresponding to the observed X-ray flux
    """
    F_nu = F_X/(nu_max-nu_min)
    # Constants
    return flx2mag(F_nu)


def EPexpo2mAB(expo,instrument:str):
    """
    expo[float, array]:     exposure time in second
    instrument[str]:        WXT or FXT
    """
    if instrument == 'WXT':
        F_X = 6.13670736e-09 * expo**-0.833761923
        E_min, E_max = 0.5, 4
        nu_max, nu_min = keV2Hz(E_max), keV2Hz(E_min)
    elif instrument == 'FXT':
        F_X = 5.49101782e-11 * expo**-0.812917739
        E_min, E_max = 0.5, 10
        nu_max, nu_min = keV2Hz(E_max), keV2Hz(E_min)
    else:
        raise KeyError('Please type correct instrument name!')
    
    return X2mAB(F_X,nu_min,nu_max)

# def register_epband():
#     'Register X-ray filter:: epwxt, epfxt'
#     wave_ = np.linspace(3.0996033431592043, 24.796826745273634,20)
#     trans_ = np.ones(20,)
#     band = sncosmo.Bandpass(wave_, trans_, name='epwxt')
#     sncosmo.register(band, 'epwxt',force=True)

#     wave_ = np.linspace(1.2398413372636818, 24.796826745273634,20)
#     trans_ = np.ones(20,)
#     band = sncosmo.Bandpass(wave_, trans_, name='epfxt')
#     sncosmo.register(band, 'epfxt',force=True)
    

def BNS_ejecta_mass(Mc, q, R1, R2, Mtov):
    pass

def NSBH_ejecta_mass(M_BH, M_NS, Chi, R_NS):
    pass

def read_curve(src,bkg,binsize=10,scale=1./12):
    with fits.open(src) as hdu:
        TSTART = hdu[0].header['TSTART']
        DATE_OBS = hdu[0].header['DATE-OBS']
        data = hdu[1].data
        TIME = data['TIME']
        RATE = data['RATE']
        ERROR = data['ERROR']
    with fits.open(bkg) as hdu:
        bkg = hdu[1].data
        TIME_bkg = bkg['TIME']
        RATE_bkg = bkg['RATE']
        ERROR_bkg = bkg['ERROR']

    t,rate,error = [],[],[]
    t_bkg,rate_bkg,error_bkg = [],[],[] #Scaled
    #Rebin
    pin = TIME[0]
    while pin+binsize < TIME[-1]:
        idx = np.where((TIME>pin) & (TIME<pin+binsize))[0]
        true_size = len(idx)
        if true_size == 0:
            pin += binsize
            continue
        else:
            t.append(sum(TIME[idx])/true_size)
            rate.append(sum(RATE[idx])/true_size)
            error.append(np.sqrt(sum(ERROR[idx]**2))/true_size)
            rate_bkg.append(scale*sum(RATE_bkg[idx])/true_size)
            error_bkg.append(scale*np.sqrt(sum(ERROR_bkg[idx]**2))/true_size)
            pin += binsize 
    return (t,rate,error), (t_bkg,rate_bkg,error_bkg), (t,np.array(rate)-np.array(rate_bkg),np.sqrt(np.array(error)**2+np.array(error_bkg)**2)), TSTART

def Txx(times,rates,c=0.9):
    """
    c[float]:   confident level [0,1], default=0.9 represent 90% credible level
    """
    accumulated_rates = accumulated_rates
    rates = rates
    rates_err = rates_err
    C5 = np.quantile(accumulated_rates,c/2)
    C95 = np.quantile(accumulated_rates,1-c/2)
    i5, i95 = np.abs(accumulated_rates-C5).argmin(), np.abs(accumulated_rates-C95).argmin()
    T5, T95 = times[i5], times[i95]
    return T95-T5

def TA_quick(obsid,snum,root='',binsize=10,pha_file=None,rebin=2,grp=False,nH=None,group=None,rx=None,sep=True,ins='WXT',plotstyle='step',N=500,get_unabs=True,chatter=10):
    """
    Perform quick analysis for TA.

    Args:
    - obsid : str : e.g. ep06800000356wxt45
    - snum : int : Source number
    - root (str, optional): Root directory of Grace data. Defaults to './'.
    - binsize (int, optional): Size of bin in seconds. Defaults to 10.
    - rebin (int, optional): Rebin factor. Defaults to 2.
    - rx (float, optional): refine x range.
    - ins (str, optional): Instrument name. Defaults to 'WXT'.
    - snum_prefix: specify snum in file searching
    """
    #Designed for FXT
    #Plot curve
    #Move to root dir
    plotmode = 'euf resid'
    os.chdir(root)
    root = ''
    snum = str(snum)
    lc_src = os.path.join(root,obsid)+snum+'.lc'
    #Designed for WXT
    #Plot curve
    arf = rmf = None
    if ins == 'WXT':
        try:
            snum = str(snum)
            lc_src = glob.glob(os.path.join(root,'**%s.lc'%snum))[0]
            lc_bkg = glob.glob(os.path.join(root,'**%sbk.lc'%snum))[0]
            pha_src = glob.glob(os.path.join(root,'**%s.pha'%snum))[0]
            pha_bkg = glob.glob(os.path.join(root,'**%sbk.pha'%snum))[0]
            arf = glob.glob(os.path.join(root,'**.arf'))[0]
            rmf = glob.glob(os.path.join(root,'**.rmf'))[0]
        except:
            print('Cannot find lc.')
    elif ins == 'FXT':
        lc_src = glob.glob(os.path.join(root,'**.lc'))[0]
        lc_bkg = glob.glob(os.path.join(root,'**bk**.lc'))[0]
        pha_src = glob.glob(os.path.join(root,'**src**.pha'))[0]
        pha_bkg = glob.glob(os.path.join(root,'**bkg**.pha'))[0]
        arf = glob.glob(os.path.join(root,'**.arf'))[0]
        rmf = glob.glob(os.path.join(root,'**.rmf'))[0]
    
    # lc_src = os.path.join(root,obsid)+snum+'.lc'
    # lc_bkg = os.path.join(root,obsid)+snum+'bk.lc'
    # pha_src = os.path.join(root,obsid)+snum+'.pha'
    # arf = os.path.join(root,obsid)+snum+'.arf'
    # rmf = os.path.join(root,obsid)+'.rmf'
    
    if pha_file:
        pha_src = pha_file
    
    erange = {'WXT':[0.5,4],'FXT':[0.5,10]}
    
    if grp and not pha_file:
        pha_grp_src = 'PC.pi'
        if os.path.exists(pha_grp_src):
            os.remove(pha_grp_src)
        print(f"grppha infile={pha_src} outfile=PC.pi chatter=0 comm='group min {group} & chkey RESPFILE {rmf} & chkey ANCRFILE {arf} & chkey BACKFILE {pha_bkg} & exit'")
        #grp_data(pha_src,outputname=pha_grp_src,arf=arf,rmf=rmf,group=group)
        os.system(f"grppha infile={pha_src} outfile=PC.pi chatter=0 comm='group min {group} & chkey RESPFILE {rmf} & chkey ANCRFILE {arf} & chkey BACKFILE {pha_bkg} & exit'") 
        pha_src = pha_grp_src

    
    models = ['powerlaw','bbody','apec']
    key_pars = {'powerlaw':['PhoIndex','lg10Flux'],'bbody':['kT','lg10Flux'],'apec':['kT','lg10Flux']}
    for model in models:
        if get_unabs:
            mname = 'tbabs*cflux*' + model
        else:
            mname = 'cflux*tbabs*' + model
        fix_ = {model+'.norm':1,'cflux.Emin':erange[ins][0],'cflux.Emax':erange[ins][1]}
        if nH:
            fix_['TBabs.nH'] = nH
        fitted_data, par_table_i = xspec_fitting(pha_src,mname=mname,grp=grp,arf=arf,rmf=rmf,rebin=rebin,instrument=ins,plotmode='ldata del',N=N,chatter=chatter,**fix_)
        
        model_leg = ''
        for key_par in key_pars[model]:
            p = par_table_i[key_par]
            p_high = par_table_i['%s_err_high'%key_par]
            p_low = par_table_i['%s_err_low'%key_par]
            
            if key_par == 'lg10Flux':
                ofm = int(np.floor(p))
                p, p_high, p_low = 10**p, 10**p_high, 10**p_low
                if get_unabs:
                    key_par = 'Unabs Flux'
                else:
                     key_par = 'Obs Flux'
                model_leg += '%s = $%.3f^{+%.3f}_{-%.3f}\\times 10^{%s}$\n'%(key_par,p/10**ofm,(p_high-p)/10**ofm,(p-p_low)/10**ofm,ofm)
            else:
                model_leg += '%s = %.3f (+%.3f/-%.3f)\n'%(key_par,p,p_high-p,p-p_low)
        if not nH:
            key_par = 'nH'
            p = par_table_i[key_par]
            p_high = par_table_i['%s_err_high'%key_par]
            p_low = par_table_i['%s_err_low'%key_par]
            model_leg += '%s = %.3e (+%.3e/-%.3e)\n'%(key_par,p,p_high-p,p-p_low)
        
        key_par = 'cstat/dof'
        p = par_table_i[key_par]
        model_leg += '%s = %.3f'%(key_par,p)
            
        xspec_plot(fitted_data,save_dir=os.path.join(root,obsid)+snum+'_%s.pdf'%model,plotstyle=plotstyle,model_leg=model_leg,title=model)
        
        par_table_i.write(os.path.join(root,'%s.csv'%model),format='csv',overwrite=True)
    
    lcurve_plot(src=lc_src,bkg=lc_bkg,binsize=binsize,save_dir=os.path.join(root,obsid)+snum+'_lc.pdf',sep=sep,rx=rx)
    
        

    
    
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
