'Cross-match GW-AGN Example'
import EPTools
from astropy.coordinates import SkyCoord
from EPTools.utils import *
import numpy as np
from astropy.cosmology import Planck18
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier, VizierClass
from ligo.skymap.postprocess import crossmatch
from ligo.skymap.io.fits import read_sky_map
from astropy.table import Table, join, join_skycoord

def crossmatch_GW_with_cat(skymap_dir,save_dir):
    """
    skymap_dir [str]:   skymap file dir or url
    save_dir[str]:      saved dir of 
    """

    #Load Cat
    wise_agn_table = Table.read("/Users/liangrunduo/EP/Catalogue/WISE_AGN.csv", format="csv")
    milliquas_table = Table.read("/Users/liangrunduo/EP/Catalogue/Milliquas.csv", format="csv")
    skymap = read_sky_map(skymap_dir,moc=True)
    #skymap = read_sky_map('/Users/liangrunduo/EP/GW/S241102br_skymap.fits',moc=True)
    milliquas_table_valid = milliquas_table[milliquas_table['Z']>0]
    dist = Planck18.luminosity_distance(milliquas_table_valid['Z'])
    coordinates = SkyCoord(milliquas_table_valid['RA']*u.deg, milliquas_table_valid['DEC']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(milliquas_table_valid[result.searched_prob_vol < 0.9])
    matched_milliquas = milliquas_table_valid[result.searched_prob_vol < 0.96]

    wise_agn_table_valid = wise_agn_table[wise_agn_table['z']>0]
    dist = Planck18.luminosity_distance(wise_agn_table_valid['z'])
    coordinates = SkyCoord(wise_agn_table_valid['_RAJ2000']*u.deg, wise_agn_table_valid['_DEJ2000']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(wise_agn_table_valid[result.searched_prob_vol < 0.9])
    matched_wise = wise_agn_table_valid[result.searched_prob_vol < 0.96]
    matched_wise.rename_column('HMQ','NAME')
    for i in range(len(matched_wise)):
        if type(matched_wise['NAME'][i]) is not np.str_:
            matched_wise['NAME'][i] = matched_wise['WISEA'][i]


    # matched_milliquas.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_milliquas.csv',format='csv',overwrite=True)
    # matched_wise.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_wise.csv',format='csv',overwrite=True)

    matched_all = join(matched_milliquas, matched_wise, keys='NAME',join_type='outer')
    print(matched_all)
    
    for i in range(len(matched_all)):
        if not isinstance(matched_all['RA'][i], np.floating):
            matched_all['RA'][i] = matched_all['_RAJ2000'][i]
            matched_all['DEC'][i] = matched_all['_DEJ2000'][i]

        if not isinstance(matched_all['Z'][i], np.floating):
            matched_all['Z'][i] = matched_all['z'][i]

    matched_all.write(save_dir,format='csv',overwrite=True)
    #print(matched_all)



crossmatch_GW_with_cat('https://gracedb.ligo.org/api/superevents/S240413p/files/bayestar.multiorder.fits',
                       '/Users/liangrunduo/EP/GW/crossmatch/S240413p_AGN.csv')