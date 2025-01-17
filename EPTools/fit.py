from .utils import *
import dynesty

def xspec_fitting(sname:str,mname:str,rebin=2,instrument:str='WXT',**fixed_par):
    if instrument == 'WXT':
        el, eh = 0.5, 4
    elif instrument == 'FXT':
        el, eh = 0.5, 10
    else:
        raise KeyError('Input Valid Instrument Type')
    
    print(fixed_par)
    
    s = Spectrum(sname)
    

def dynesty_fitting():
    pass