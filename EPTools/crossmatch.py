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


    