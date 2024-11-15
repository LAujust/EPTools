import EPTools
import EPTools.gw_em
import EPTools.utils
import numpy as np
import pandas as pd


if __name__ == '__main__':

    crossmatch = EPTools.Crossmatch()
    file_root = '/Users/liangrunduo/EP/GW/Plans/S241109bn/'
    wxt_fname = 'wxt_obs_S241109bn.csv'
    fxt_fname = 'fxt_obs_S241109bn.csv'



    EPTools.gw_em.make_WXTPlan(file_root+wxt_fname)


    # asassn_lc = crossmatch.asassn_lc(pos='18:54:11.5 -88:02:55.22',radius=50,units='arcsec')
    # print(asassn_lc)
    # print(asassn_lc.data)

    #pos = '28.8473 5.9364'
    #print(crossmatch.ztf_lc(pos,band='r',match_rad=5))
    #table.pprint()