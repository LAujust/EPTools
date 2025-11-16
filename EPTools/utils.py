import numpy as np
import pandas as pd
import sys, os, glob, re
import astropy.units as u
import subprocess, warnings
import astropy.constants as c
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table
from astropy.utils.data import conf
from ligo.gracedb.rest import GraceDb
import ligo.skymap
from gdpyc import GasMap, DustMap
from scipy.special import gammainc, gammaincc, gammaincinv
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.time import Time
warnings.filterwarnings("ignore")
sys.path.append('/Users/liangrunduo/heasoft-6.34/aarch64-apple-darwin23.5.0/lib/python')
sys.path.append('$HEADAS/lib/python')
import xspec as xs
try:
    import xspec
except ImportError:
    print("HEASoft is not initialized or heasoftpy not installed. ")
    class _HeasoftDummy:
        def __getattr__(self, name):
            raise RuntimeError(
                "HEASoft is not initialized or heasoftpy not installed. "
                "Please run 'heasoft' and try again."
            )
    heasoftpy = _HeasoftDummy()
    
from .plot import *
from .fit import *




class HeaEnv:
    def __init__(self,headas=None,caldb=None):
        self.headas = headas
        self.caldb = None
    
    def init_in_shell():
        pass
    
    def init_in_notebook():
        pass
    
    def _search_init_in_file():
        pass


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

def flx2lum(f:float,d:float):
    """convert observed flux to luminosity

    Args:
        f (float): flux in [erg/cm^2/s]
        d (float): distance in [pc]

    Returns:
        float: luminosity in [erg/s]
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

def li_ma_sigma(N_on:float, N_off:float, alpha:float):
    """significance based on Li-Ma formula

    Args:
        N_on (float): photons within src region
        N_off (float): photons within bkg region
        alpha (float): area ratio

    Returns:
        float: li-ma significance
    """
    if N_on <= 0 or N_off <= 0:
        return 0.0  # or np.nan
    term1 = N_on * np.log((1 + alpha) / alpha * (N_on / (N_on + N_off)))
    term2 = N_off * np.log((1 + alpha) * (N_off / (N_on + N_off)))
    return np.sqrt(2 * (term1 + term2))

def X_UL(Nsrc,Nbkg,exposure,alpha=1/12,factor=1e-9,CL = 0.9):
    """X-ray upper limit

    Args:
        Nsrc (float): photons in src region
        Nbkg (float): photons in bkg region
        exposure (float): exposure time in [s]
        alpha (float, optional): area ratio. Defaults to 1/12.
        factor (float, optional): ctr to flux convertion factor. Defaults to 1e-9.
        CL (float, optional): credible level. Defaults to 0.9.

    Returns:
        float: upper limit in CL
    """
    
    B = Nbkg * alpha
    C_1 = gammaincc(Nsrc+1,B)
    sC_1 = gammainc(Nsrc+1,B)
    part1 = gammaincinv(Nsrc+1,CL*C_1+sC_1)
    UL = part1 - B
    return UL * factor / exposure
    
def get_ctrt_to_flux(source_spec, energy_l, energy_h, nH_Galactic, PhoIndex, get_unabs=True,nH_Intrinsic=0, ins='WXT'):
    """estimate ctr to flux convertion factor
    Args:
        source_spec (str): source spectrum filename
        energy_l (float): enegry lower bound
        energy_h (float): enegry higher bound
        nH (float): nH in tbabs model
        PhoIndex (float): phoindex of powerlaw model
        get_unabs (bool, optional): calculate absorbed flux. Defaults to True.

    Returns:
        float: convertion factor
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

def read_curve(src:str,bkg:str,binsize:int=10,scale:float=1./12):
    """read .lc file

    Args:
        src (str): src curve
        bkg (str): bkg curve
        binsize (int, optional): time bin size. Defaults to 10.
        scale (float, optional): area ratio. Defaults to 1/12.

    Returns:
        tuple: src, bkg rate and TSTART
    """
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


def check_cosmic_ray(path:str,s_num:int,cts_thres:float=5):
    """check whether the detection is cosmic ray

    Args:
        path (str): path to your data (standard WXT level 2 data)
        s_num (int): source number
        cts_thres (float, optional): count rate threshold. Defaults to 5.
    """
    
    cat_hdu = fits.open(glob.glob(os.path.join(path,'**.cat')))
    evt_hdu = fits.open(glob.glob(os.path.join(path,'**po_cl.evt')))
    uf_hdu = fits.open(glob.glob(os.path.join(path,'**po_uf.evt')))

    hdudata = cat_hdu[1].data
    x0 = hdudata['x'][s_num-1]
    y0 = hdudata['y'][s_num-1]

    x0 = x0 * 16
    y0 = y0 * 16
    hdudata = evt_hdu[1].data
    data_cl = np.transpose(np.vstack((hdudata['time'], hdudata['rawx'], hdudata['rawy'], hdudata['x'], hdudata['y'], hdudata['pi'])))
    data_cl = pd.DataFrame(data_cl, columns=['time', 'rawx', 'rawy', 'x', 'y', 'pi'])

    hdudata = uf_hdu[1].data
    data_uf = np.transpose(np.vstack((hdudata['time'], hdudata['x'], hdudata['y'], hdudata['pi'], hdudata['cmosfram'])))
    data_uf = pd.DataFrame(data_uf, columns=['time', 'x', 'y', 'pi', 'cmosfram'])

    data = pd.merge(data_cl, data_uf, how='left', on=['time', 'x', 'y', 'pi'])

    fig1 = plt.figure()
    plt.hist2d(data['x'], data['y'], bins=[200, 200])     # 探测图像，可不看，仅供参考

    frame = data['cmosfram'].values
    x = data['x'].values
    y = data['y'].values
    fig2 = plt.figure()
    frame_index = np.arange(np.min(frame), np.max(frame)+1, 1)
    data_hist = np.histogram(frame, bins=frame_index)[0]
    data_hist_all = data_hist
    plt.plot(frame_index[:-1]-frame_index[0], data_hist, color='blue')
#    plt.xlim(0,100000)
    plt.ylim(0,20)
    index = (np.abs(x-x0)<50) & (np.abs(y-y0)<50)
    rawx0 = np.mean(data['rawx'].values[index])
    rawy0 = np.mean(data['rawy'].values[index])
    data_hist = np.histogram(frame[index], bins=frame_index)[0]
    plt.plot(frame_index[:-1]-frame_index[0], data_hist, color='red')
    plt.xlabel('Frame')
    plt.ylabel('Counts')
    # 全帧以及源的帧计数光变曲线，判断源的光变是否有单帧粒子数的增加引起

    # 提取单帧图像，阈值这里设为20，具体参考上一步源光变曲线的计数
    index_all = np.where(data_hist > cts_thres)[0]
    if len(index_all)>0:
        for i in index_all:
            index = (frame == frame_index[i])

            fig3 = plt.figure()
            rawx = data['rawx'][index]
            rawy = data['rawy'][index]
            plt.hist2d(rawx, rawy, bins=[50, 50], range=[[rawx0-200, rawx0+200], [rawy0-200, rawy0+200]])  # 源图像
            fig4 = plt.figure()
            plt.hist2d(rawx, rawy, bins=[100, 100], range=[[0, 4096], [0, 4096]])  # 全帧图像

    plt.show()


def TA_quick(obsid:str,snum:int,root:str='',binsize:int=10,pha_file=None,rebin=2,grp:bool=False,nH=None,group=None,rx=None,sep:bool=True,ins:str='WXT',plotstyle:str='step',N:int=500,get_unabs=True,chatter:int=10,module:str='B'):
    """Quick assembled tool for temporal and spectrum analysis for TA.

    Args:
        obsid (str): Observation ID;
        snum (int): source number;
        root (str, optional): path to data folder. Defaults to '';
        binsize (int, optional): lightcurve bin size. Defaults to 10.
        pha_file (optional): pha file name if you already have a grouped pha file. Defaults to None;
        rebin (int, optional): xspec plot rebin. Defaults to 2;
        grp (bool, optional): whether to group data. Defaults to False;
        nH (_type_, optional): Galactic hydrogen column density, take as free parameter if set None. Defaults to None;
        group (optional): minimal photons in a group, passing to grppha. Defaults to None;
        rx (optional): refined x range in curve plot. Defaults to None;
        sep (bool, optional): whether to plot src/net/bkg curve in different subplots. Defaults to True;
        ins (str, optional): EP instruments data to be proceeded. Defaults to 'WXT';
        plotstyle (str, optional): plot spectrum style of model. Defaults to 'step';
        N (int, optional): Number of iteration for fitting. Defaults to 500;
        get_unabs (bool, optional): calculate unabsorbed flux. Defaults to True;
        chatter (int, optional): chatter level. Defaults to 10.
        module (str, optional): FXT module (A or B). Defaults to 'B'.

    Raises:
        KeyError: _description_
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
        if module == None:
            try:
                pha_src = glob.glob(os.path.join(root,'fxt_b**src**.pha'))[0]
                pha_bkg = glob.glob(os.path.join(root,'fxt_b**bkg**.pha'))[0]
                arf = glob.glob(os.path.join(root,'fxt_b**src**.arf'))[0]
                rmf = glob.glob(os.path.join(root,'fxt_b**src**.rmf'))[0]
            except:
                try:
                    pha_src = glob.glob(os.path.join(root,'fxt_a**src**.pha'))[0]
                    pha_bkg = glob.glob(os.path.join(root,'fxt_a**bkg**.pha'))[0]
                    arf = glob.glob(os.path.join(root,'fxt_a**src**.arf'))[0]
                    rmf = glob.glob(os.path.join(root,'fxt_a**src**.rmf'))[0]
                except:
                    pha_src = glob.glob(os.path.join(root,'**src**.pha'))[0]
                    pha_bkg = glob.glob(os.path.join(root,'**bkg**.pha'))[0]
                    arf = glob.glob(os.path.join(root,'**.arf'))[0]
                    rmf = glob.glob(os.path.join(root,'**.rmf'))[0]
        else:
            try:
                pha_src = glob.glob(os.path.join(root,'fxt_%s**src**.pha'%(module.lower())))[0]
                pha_bkg = glob.glob(os.path.join(root,'fxt_%s**bkg**.pha'%(module.lower())))[0]
                arf = glob.glob(os.path.join(root,'fxt_%s**src**.arf'%(module.lower())))[0]
                rmf = glob.glob(os.path.join(root,'fxt_%s**src**.rmf'%(module.lower())))[0]
            except:
                raise KeyError('No valid Files!')
        
        try:
            lc_src = glob.glob(os.path.join(root,'**.lc'))[0]
            lc_bkg = glob.glob(os.path.join(root,'**bk**.lc'))[0]
        except Exception as e:
            print('Cannot find lc: %s'%e)
            lc_src = lc_bkg = None
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
            
        xspec_plot(fitted_data,save_dir=os.path.join(root,obsid)+snum+'_%s.png'%model,plotstyle=plotstyle,model_leg=model_leg,title=model)
        
        par_table_i.write(os.path.join(root,'%s.csv'%model),format='csv',overwrite=True)
    
    lcurve_plot(src=lc_src,bkg=lc_bkg,binsize=binsize,save_dir=os.path.join(root,obsid)+snum+'_lc.png',sep=sep,rx=rx)
    
        

    
    
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
