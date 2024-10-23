import EPTools
import EPTools.utils
import numpy as np
import pandas as pd



if __name__ == '__main__':
    f = 3.3e-10 
    d = 100 * 1e6
    print(EPTools.utils.flx2lum(f,d))

    crossmatch = EPTools.Crossmatch()
    pos = '12h29m06.70s +02d03m08.7s'
    data = crossmatch.xmm_archive(pos=pos,radius='1 degree')
    print(data[0]['FLUX_B8'])
    print(data[0]['FLUX_B8_ERROR'])

    asassn_lc = crossmatch.asassn_lc(pos='18:54:11.5 -88:02:55.22',radius=50,units='arcsec')
    print(asassn_lc)
    print(asassn_lc.data)
    #table.pprint()