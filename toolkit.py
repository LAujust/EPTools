#aujust@mail.ustc.edu.cn
#Toolkit for Aujust

import numpy as np
import sncosmo
import astropy.units as u
import astropy.constants as c
import os
import re
import pandas as pd
import requests
from bs4 import BeautifulSoup
from astropy.time import Time
from ligo.skymap.io.fits import read_sky_map
import healpy as hp
import pickle
import scipy

def __check__():
    print('Check~')


def load_filter_list():
    all = True
    filt_list = ['sdssu','sdssg','sdssr','sdssi','sdssz','desg','desr','desi','desz','desy','f435w','f475w','f555w','f606w','f625w',
        'f775w','nicf110w','nicf160w','f098m','f105w','f110w','f125w','f127m','f139m','f140w','f153m','f160w','f218w','f225w',
        'f275w','f300x','f336w','f350lp','f390w','f689m','f763m','f845m','f438w','uvf555w','uvf475w',
        'uvf606w','uvf625w','uvf775w','uvf814w','uvf850lp','cspb','csphs','csphd','cspjs','cspjd','cspv3009',
        'cspv3014','cspv9844','cspys','cspyd','cspg','cspi','cspk','cspr','cspu','f070w','f090w','f115w','f150w',
        'f200w','f277w','f356w','f444w','f140m','f162m','f182m','f210m','f250m','f300m','f335m','f360m','f410m','f430m',
        'f460m','f480m','lsstu','lsstg','lsstr','lssti','lsstz','lssty','keplercam::us','keplercam::b','keplercam::v','keplercam::v',
        'keplercam::r','keplercam::i','4shooter2::us','4shooter2::b','4shooter2::v','4shooter2::r','4shooter2::i','f062','f087',
        'f106','f129','f158','f184','f213','f146','ztfg','ztfr','ztfi','uvot::b','uvot::u','uvot::uvm2','uvot::uvw1','uvot::uvw2',
        'uvot::v','uvot::white','ps1::open','ps1::g','ps1::r','ps1::i','ps1::z','ps1::y','ps1::w','atlasc','atlaso','2massJ',
        '2massH','2massKs','wfst_u','wfst_g','wfst_r','wfst_i','wfst_z','wfst_w','bg_q'
    ]
    for filt in filt_list:
        try:
            _x = sncosmo.get_bandpass(filt)
        except:
            print('Fail for '+filt)
            if all:
                all = False
    if all:
        print('Load all filters successfully!')
    return filt_list

def load_wfst_bands():
    add_bands = ['u','g','r','i','w','z']
    wfst_bands = ['wfst_'+i for i in add_bands]
    try:
        for add_band in add_bands:
            data = np.loadtxt('/home/Aujust/data/Kilonova/WFST/transmission/WFST_WFST.'+add_band+'_AB.dat')
            wavelength = data[:,0]
            trans = data[:,1]
            band = sncosmo.Bandpass(wavelength, trans, name='wfst_'+add_band)
            sncosmo.register(band, 'wfst_'+add_band)
            
        #Load BlackGEM/MeerLICHT filter
        data = np.loadtxt('/home/Aujust/data/Kilonova/WFST/transmission/BlackGEM/BlackGEM_q.dat')
        wavelength = data[:,0]
        trans = data[:,1]
        band = sncosmo.Bandpass(wavelength, trans, name='bg_q')
        sncosmo.register(band, 'bg_q')
        print('Sucess.')
        
    except:
        pass
    return wfst_bands

def mab2flux(mab):
    #erg s^-1 cm^-2
    return 10**(-(mab+48.6)/2.5)

def flux2mab(f):
    #erg s^-1 cm^-2
    return -2.5*np.log10(f)-48.6

def sumab(mab_list):
    _flux_all = 0
    for mab in mab_list:
        _flux_all += mab2flux(mab)
    return flux2mab(_flux_all)

def get_all_event_name(event_folder_path):
    pattern = re.compile(r'S[^.]+')
    event_names = []
    filenames = os.listdir(event_folder_path)
    for filename in filenames:
        event_names.append(pattern.search(filename).group(0)[:])
    return event_names

def event_mjd(event_id):
    'event_id[str]: ID should be conicide with IDs in GraceDb and return the trigger time in MJD'
    url = 'https://gracedb.ligo.org/superevents/{}/view/#event-information'.format(event_id)
    response = requests.get(url)
    soup = BeautifulSoup(response.text,'html.parser')
    html_str = str(soup.find_all('time')[1])
    soupt = BeautifulSoup(html_str,'html.parser')
    return Time(soupt.time['utc'][:-4],format='iso',scale='utc').mjd

def load_skymap(skymap_dir,multiorder=False):
    if multiorder:
        m, meta = read_sky_map(skymap_dir,nest=False,distances=True)
        prob, distmu, distsigma = m[0], m[1], m[2]
    else:
        prob, distmu,distsigma,distnorm = hp.fitsfunc.read_map(skymap_dir,field=[0,1,2,3])
        
    skymap = {
                'prob':prob,
                'distmu':distmu,
                'distsigma':distsigma
                }
    return skymap

def prior_mtov_Tim(m):
    '''
    Construct Likelihood/Prob of Mtov reffering Tim 2020.
    '''
    mean_sig = np.array([
    [2.14,0.1],
    [2.01,0.04],
    [1.908,0.016],
    [2.16,0.17]
    #[2.13,0.09]
    ])

    li = 1
    for i in range(mean_sig.shape[0]):
        x = (m - mean_sig[i,0])/mean_sig[i,1]
        li_i = scipy.stats.norm.cdf(x)
        if i == 3:
            li = li*(1-li_i)
        else:
            li = li*li_i

    return li

def blackbody(params,lam):
    T = params[0]  #K
    R = params[1]  #cm
    #print(T,R,lam)
    const = 8*pi**2*R**2*H_CGS*C_CGS**2*ANG_CGS/lam**5
    const2 = H_CGS*C_CGS/K_B_CGS
    return const/(np.exp(const2/(lam*T))-1)  #erg s^-1 A-1

def wfst_obs_to_plan(wfst_obs,save=False,output_name='wfst_plan'):
    import simsurvey
    
    wfst_survey = pd.read_csv(wfst_obs)

    'drop invalid data'
    fail_index = wfst_survey[np.isnan(wfst_survey['maglim'])]
    fail_index = fail_index.index
    wfst_survey = wfst_survey.drop(fail_index)
    
    try: 
        ras = wfst_survey['ra']
        decs = wfst_survey['dec']
    except:
        from astropy.coordinates import SkyCoord
        ras = []
        decs = []
        for coord in wfst_survey['Ra,Dec']:
            ra,dec = coord.split()
            cod = SkyCoord(ra,dec,frame='icrs')
            ras.append(cod.ra.degree)
            decs.append(cod.dec.degree)

    wfst_survey['ra'] = ras
    wfst_survey['dec'] = decs
    wfst_survey['band'] = ['wfst_'+wfst_survey.loc[i,'Filter'] for i in list(wfst_survey.index)]
    wfst_survey['time'] = wfst_survey['mjd']

    obs = {'time': [], 'ra': [], 'dec': [], 'band': [], 'maglim': [], 'skynoise': [], 'comment': [], 'zp': []}

    for k in wfst_survey.keys():
        obs[k] = wfst_survey[k]

    obs['zp'] = [26 for i in list(wfst_survey.index)]
    obs['comment'] = ['' for i in list(wfst_survey.index)]        
    obs['skynoise'] = 10**(-0.4 * (np.array(obs['maglim']) - 26)) / 5

    'Create simsurvey.SurveyPlan object'
    plan = simsurvey.SurveyPlan(time=obs['time'],
                                    band=obs['band'],
                                    ra=obs['ra'],
                                    dec=obs['dec'],
                                    skynoise=obs['skynoise'],
                                    obs_ccd=None,
                                    zp=obs['zp'],
                                    comment=obs['comment'],
                                    height=2.6,
                                    width=2.6
                                    )

    print(plan.pointings)

    if save:
        with open(output_name+'.pkl','wb') as f:
            pickle.dump(plan,f)
        f.close()
        
    return plan

def load_skymap(skymap_dir,multiorder=False):
    if multiorder:
        m, meta = read_sky_map(skymap_dir,nest=False,distances=True)
        prob, distmu, distsigma = m[0], m[1], m[2]
    else:
        prob, distmu,distsigma,distnorm = hp.fitsfunc.read_map(skymap_dir,field=[0,1,2,3])
        
    skymap = {
                'prob':prob,
                'distmu':distmu,
                'distsigma':distsigma
                }
    return skymap

#--------------------------------------#
#               Constant               #
#--------------------------------------#

day2sec = u.day.cgs.scale
MPC_CGS = u.Mpc.cgs.scale
C_CGS = c.c.cgs.value
M_SUN_CGS = c.M_sun.cgs.value
G_CGS = c.G.cgs.value
Jy = u.Jy.cgs.scale
ANG_CGS = u.Angstrom.cgs.scale
pi = np.pi
pc10 = 10 * u.pc.cgs.scale
SB_CGS = c.sigma_sb.cgs.value
H_CGS = c.h.cgs.value
K_B_CGS = c.k_B.cgs.value
KM_CGS = u.km.cgs.scale