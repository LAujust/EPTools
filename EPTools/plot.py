from .utils import *


def plot_gcn_data(file_root='./',output='standard_gcn.pdf'):
    '''
    Plot GCN Circilar data and EP-WXT/FXT data from .csv file
    
    Data should be stored in .csv file 

    Data Format:
    
    x_data.csv: [dt, instrument, flux, fluxerr, detection]
    optical_data: [dt, telescope, mag, magerr, band, detection]
    radio_data: [dt, instrument, flux, fluxerr, detection, frequency]  
    #flux in radio data often refers to flux density [mJy/uJy]
    '''
    
    csv_files = glob.glob("*.csv")

