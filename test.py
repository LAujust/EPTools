import EPTools
import EPTools.utils
import numpy as np
import pandas as pd
from astropy.utils.data import conf
from astropy.config import create_config_file
from astropy.config import get_config_dir
#create_config_file('astropy') 
print(get_config_dir())   #You can create a config file and modify 'timeout' parameters to avoid readout timeout. 


if __name__ == '__main__':

    f = 3.3e-10 
    d = 100 * 1e6
    print(EPTools.utils.flx2lum(f,d))

    crossmatch = EPTools.Crossmatch()
    pos = '12h29m06.70s +02d03m08.7s'
    data = crossmatch.xmm_archive(pos=pos,radius='1 degree')
    print(data[0]['FLUX_B8'])
    print(data[0]['FLUX_B8_ERROR'])

    # asassn_lc = crossmatch.asassn_lc(pos='18:54:11.5 -88:02:55.22',radius=50,units='arcsec')
    # print(asassn_lc)
    # print(asassn_lc.data)

    ra, dec = 28.8473, 5.9364
    print(crossmatch.ztf_lc(ra=ra, dec=dec,band='r',match_rad=5))
    #table.pprint()