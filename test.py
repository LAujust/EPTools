import EPTools
import EPTools.utils
import numpy as np
import pandas as pd


f = 3.3e-10 
d = 100 * 1e6
print(EPTools.utils.flx2lum(f,d))

crossmatch = EPTools.Crossmatch()
pos = '12h29m06.70s +02d03m08.7s'
data = crossmatch.XMM_archive(pos=pos,radius='1 degree')
print(data[0]['FLUX_B8'])
print(data[0]['FLUX_B8_ERROR'])

#table.pprint()