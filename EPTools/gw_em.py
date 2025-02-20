'GW-EM'
from .utils import *
sys.path.append('$HEADAS/lib/python')
from xspec import * 
'A module to transform synthetic model SED file to observed SED of X-ray detector'
'With the help of PyXspec'

'Register X-ray filter:: epwxt, epfxt'
def register_epband():
    wave_ = np.linspace(3.0996033431592043, 24.796826745273634,20)
    trans_ = np.ones(20,)
    band = sncosmo.Bandpass(wave_, trans_, name='epwxt')
    sncosmo.register(band, 'epwxt',force=True)

    wave_ = np.linspace(1.2398413372636818, 24.796826745273634,20)
    trans_ = np.ones(20,)
    band = sncosmo.Bandpass(wave_, trans_, name='epfxt')
    sncosmo.register(band, 'epfxt',force=True)

def observed_SED():
    AllData.clear()
    pass

def make_EPPlan(filename:str,instrument:str):
    zp = 25.
    pointings = pd.read_csv(filename)
        
    if instrument == 'WXT':
        width, height = np.sqrt(75), np.sqrt(75)
        band = 'epwxt'
        ra = pointings['Pointing RA']
        dec = pointings['Pointing Dec']
        expo = pointings['Exposure (s)']
        time = Time(list(pointings['Obs Start Time (UTC)'].to_numpy()),format='isot', scale='utc').mjd

    elif instrument == 'FXT':
        width, height = 1, 1
        band = 'epfxt'
        ra = pointings['Obj_RA']
        dec = pointings['Obj_Dec']
        expo = pointings['Exposure Time (s)']
        time = Time(list(pointings['Obs Start Time (UTC) '].to_numpy()), scale='utc').mjd
    

    mAB = EPexpo2mAB(expo,instrument=instrument)
    zp = zp * np.ones(ra.shape)
    band = [band] * len(ra)
    skynoise = 10**(-0.4 * (mAB - zp)) / 5

    plan = simsurvey.SurveyPlan(time=time,
                                band=band,
                                ra=ra,
                                dec=dec,
                                skynoise=skynoise,
                                obs_ccd=None,
                                zp=zp,
                                comment=[' ']*len(ra),
                                height=height,
                                width=width
                                )
    
    return plan

    


