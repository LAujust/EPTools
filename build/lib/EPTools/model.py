from .utils import *


class Model(object):
    def __init__(self) -> None:
        pass

    def calc_lightcurve(self):
        pass

    def calc_spectra(self):
        pass

class BlackBody(Model):
    def __init__(self):
        pass
    def calc_spectra(self,T,nus=np.logspace(14,19,200)):
        'T[float]:  temperature in K'
        h = c.h
        C = c.c
        k = c.k_B
        nus = nus*u.Hz
        T = T*u.K
        sed = 2*h*nus**3/(C**2*(np.exp((h*nus/(k*T)).cgs)-1))
        return sed.cgs


class Sun19(Model):
    def __init__(self) -> None:
        pass