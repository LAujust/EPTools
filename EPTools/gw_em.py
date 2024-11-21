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

def make_WXTPlan(wxt_filename:str,event_id='MS'):

    wxt_pointings = pd.read_csv(wxt_filename)
    print(wxt_pointings)

    


