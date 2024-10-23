from .utils import *

class Crossmatch(object):
    def __init__(self):
        self.info = None
        self.heasarc = Heasarc()
        self.asassn = SkyPatrolClient()

    

    def xmm_archive(self,pos,radius=3*u.arcmin):
        """
        pos[str]:       astropy.SkyCoords input, i.e. 'ra dec'(in degree), 'hhmmss +ddmmss', 'hh:mm:ss +dd:mm:ss'
        radius[units]:  arcmin/arcsec/degree (i.e. u.arcmin)
        """
        coord = SkyCoord(pos, frame='icrs')
        table = self.heasarc.query_region(coord, mission='xmmslewful', radius=radius)
        return table
    
    def asassn_lc(self,pos,radius=3,units='arcmin',catalog='master_list'):
        """
        pos[str]:       'ra_deg dec_deg' or 'hh:mm:ss +dd:mm:ss'
        radius[float]:  radius in degree
        """
        string = pos.split(' ')
        ra, dec = string[0],string[1]
        lcs = self.asassn.cone_search(ra_deg=ra, dec_deg=dec, 
                                      radius=radius, units=units, download=True, 
                                      #save_dir='./', file_format='csv',
                                      catalog=catalog, threads=8)

        return lcs
    
    def ztf_lc(self, ra, dec, band:str, match_rad:float=5):
        """
        :param band: str
        :param match_rad: float, unit is arcsec; defaule value is 5 arcsec
        :return: astropy.io.votable
        """
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

    