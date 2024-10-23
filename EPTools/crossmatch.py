from .utils import *

class Crossmatch(object):
    def __init__(self):
        self.info = None
        self.heasarc = Heasarc()

    def XMM_archive(self,pos,radius=3*u.arcmin):
        coord = SkyCoord(pos, frame='icrs')
        table = self.heasarc.query_region(coord, mission='xmmslewful', radius=radius)
        return table


    