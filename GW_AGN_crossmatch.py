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
from astropy.table import Table

def crossmatch_GW_with_cat(skymap_dir):
    """
    skymap_dir [str]:   skymap file dir or url
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
    print(milliquas_table_valid[result.searched_prob_vol < 0.9])

    wise_agn_table_valid = wise_agn_table[wise_agn_table['z']>0]
    dist = Planck18.luminosity_distance(wise_agn_table_valid['z'])
    coordinates = SkyCoord(wise_agn_table_valid['_RAJ2000']*u.deg, wise_agn_table_valid['_DEJ2000']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    print(wise_agn_table_valid[result.searched_prob_vol < 0.9])


crossmatch_GW_with_cat('https://gracedb.ligo.org/api/superevents/S241102br/files/bayestar.multiorder.fits')