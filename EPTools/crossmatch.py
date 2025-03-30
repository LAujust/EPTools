from .utils import *
from astropy.cosmology import Planck18
from astroquery.vizier import Vizier, VizierClass
from ligo.skymap.postprocess import crossmatch
from ligo.skymap.io.fits import read_sky_map
from astropy.table import Table, join, join_skycoord, vstack

class Crossmatch(object):
    def __init__(self):
        self.info = None
        self.heasarc = Heasarc()
        self.asassn = None
        self.recommand_xmatch_mission = ['swiftmastr','rassmaster','xmmslewful','xmmcdfs210',
                                         'xmmcdfs510','xmmssc','xmmssclwbs','xmmstack',
                                         'xmmstackob',
                                         ]
    
    
    def xmm_slew_archive(self,pos,radius=3*u.arcmin):
        return self.xmatch_archive(pos=pos,radius=radius,mission='xmmslewful')

    def rosat_archive(self,pos,radius=3*u.arcmin):
        return self.xmatch_archive(pos=pos,radius=radius,mission='rassmaster')
    
    def swift_archive(self,pos,radius=3*u.arcmin):
        return self.xmatch_archive(pos=pos,radius=radius,mission='swiftmastr')
    
    def asassn_lc(self,pos,radius=3,units='arcmin',catalog='master_list'):
        """
        pos[str]:       'ra_deg dec_deg' or 'hh:mm:ss +dd:mm:ss'
        radius[float]:  radius in degree
        """
        if self.asassn is None:
            self.asassn = SkyPatrolClient()
        string = pos.split(' ')
        ra, dec = string[0],string[1]
        lcs = self.asassn.cone_search(ra_deg=ra, dec_deg=dec, 
                                      radius=radius, units=units, download=True, 
                                      #save_dir='./', file_format='csv',
                                      catalog=catalog, threads=8)

        return lcs
    
    def ztf_lc(self, pos:str, band:str, match_rad:float=5):
        """
        :param band: str
        :param match_rad: float, unit is arcsec; defaule value is 5 arcsec
        :return: astropy.io.votable
        """
        ra, dec = pos.split(' ')
        ra, dec = float(ra), float(dec)
        match_rad = match_rad/3600
        API =  f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{ra}%20{dec}%20{match_rad}&BANDNAME={band}"
        lcdata = parse_single_table(API).to_table()
        return lcdata
    
    def crossmatch_milliquas(self,ra,dec,radius):
        """
        ra, dec: in deg
        radius: in arcsec
        """
        #catalog_mlq = np.loadtxt('/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/milliquas/milliquas.txt')
        #catalog_mlq = ascii.read('/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/milliquas/milliquas.txt',format='txt')
        catalog_mlq = pd.read_csv('/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/milliquas/milliquas.txt',
                                on_bad_lines='skip')
        print(catalog_mlq)

    

    def xmatch_archive(self,pos,radius=3*u.arcmin,mission='xmmslewful'):
        """
        pos[str]:       astropy.SkyCoords input, i.e. 'ra dec'(in degree), 'hhmmss +ddmmss', 'hh:mm:ss +dd:mm:ss'
        radius[units]:  arcmin/arcsec/degree (i.e. u.arcmin)
        """
        coord = SkyCoord(pos, frame='icrs')
        table = self.heasarc.query_region(coord, mission='xmmslewful', radius=radius)
        return table
    
def crossmatch_GW_with_cat(skymap_dir,save_dir):
    """
    skymap_dir [str]:   skymap file dir or url
    save_dir[str]:      saved dir of 
    """

    #Load Cat
    wise_agn_table = Table.read("/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/WISE_AGN.csv", format="csv")
    milliquas_table = Table.read("/Users/liangrunduo/Desktop/Aujust/NAOC/EP/Crossmatch/catalogs/Milliquas.csv", format="csv")
    skymap = read_sky_map(skymap_dir,moc=True)
    #skymap = read_sky_map('/Users/liangrunduo/EP/GW/S241102br_skymap.fits',moc=True)
    milliquas_table_valid = milliquas_table[milliquas_table['Z']>0]
    dist = Planck18.luminosity_distance(milliquas_table_valid['Z'])
    coordinates = SkyCoord(milliquas_table_valid['RA']*u.deg, milliquas_table_valid['DEC']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(milliquas_table_valid[result.searched_prob_vol < 0.9])
    matched_milliquas = milliquas_table_valid[result.searched_prob_vol < 0.95]

    wise_agn_table_valid = wise_agn_table[wise_agn_table['z']>0]
    dist = Planck18.luminosity_distance(wise_agn_table_valid['z'])
    coordinates = SkyCoord(wise_agn_table_valid['_RAJ2000']*u.deg, wise_agn_table_valid['_DEJ2000']*u.deg, dist)
    result = crossmatch(skymap, coordinates)
    #print(wise_agn_table_valid[result.searched_prob_vol < 0.9])
    matched_wise = wise_agn_table_valid[result.searched_prob_vol < 0.95]
    matched_wise.rename_column('HMQ','NAME')
    for i in range(len(matched_wise)):
        if type(matched_wise['NAME'][i]) is not np.str_:
            matched_wise['NAME'][i] = matched_wise['WISEA'][i]


    # matched_milliquas.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_milliquas.csv',format='csv',overwrite=True)
    # matched_wise.write('/Users/liangrunduo/EP/GW/crossmatch/S240413p_wise.csv',format='csv',overwrite=True)

    if len(matched_milliquas) * len(matched_wise) > 0:
        matched_all = join(matched_milliquas, matched_wise, keys='NAME',join_type='outer')
    else:
        matched_all = vstack([matched_milliquas,matched_wise])
    print(matched_all)
    
    for i in range(len(matched_all)):
        if not isinstance(matched_all['RA'][i], np.floating):
            matched_all['RA'][i] = matched_all['_RAJ2000'][i]
            matched_all['DEC'][i] = matched_all['_DEJ2000'][i]

        if not isinstance(matched_all['Z'][i], np.floating):
            matched_all['Z'][i] = matched_all['z'][i]

    matched_all.write(save_dir,format='csv',overwrite=True)


def match_cat(source_cat,cat,radius,known_source_cat=None,nthneighbor=1):
    if known_source_cat:
        r = 3.5 * u.arcmin
        idx, sep, _ = source_cat.match_to_catalog_sky(known_source_cat,nthneighbor=nthneighbor)
        known_source_idx = sep < r
        known_source_cat_idx = idx[known_source_idx]
    else:
        known_source_idx, known_source_cat_idx = [], []

    idx, sep, _ = source_cat.match_to_catalog_sky(cat,nthneighbor=nthneighbor)
    filtered_id = sep < radius
    cat_matched_idx, cat_matched_sep = idx[filtered_id], sep[filtered_id]
    source_matched_idx = filtered_id
    return source_matched_idx, cat_matched_idx, known_source_idx, known_source_cat_idx