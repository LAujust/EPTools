'GW-EM'
from .utils import *
sys.path.append('$HEADAS/lib/python')
from xspec import *
'A module to transform synthetic model SED file to observed SED of X-ray detector'
'With the help of PyXspec'

'Register X-ray filter:: epwxt, epfxt'
def register_epband():
    wave_ = np.linspace(Hz2lam(keV2Hz(4)), Hz2lam(keV2Hz(0.5)),20)
    trans_ = np.ones(20,)
    sncosmo.Bandpass(wave_, trans_, name='epwxt')

    wave_ = np.linspace(Hz2lam(keV2Hz(10)), Hz2lam(keV2Hz(0.5)),20)
    trans_ = np.ones(20,)
    sncosmo.Bandpass(wave_, trans_, name='epfxt')

def observed_SED():
    AllData.clear()
    pass

def make_WXTPlan(wxt_filename:str,event_id='MS'):

    wxt_pointings = pd.read_csv(wxt_filename)
    print(wxt_pointings)

    

register_epband()

