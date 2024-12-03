from .utils import *


def plot_gcn_data(file_dir='',output='standard_gcn.pdf'):
    '''
    Plot GCN Circilar data and EP-WXT/FXT data from .csv file

    Data Format:
    
    data.csv: [dt, instrument, mag, magerr, flux, fluxerr, detection, type] 
    #flux in radio data often refers to flux density [mJy/uJy]
    #type should be either optical, x-ray or radio
    '''

    Alldata = pd.read_csv(file_dir).sort_values(by='dt').reset_index(drop=True)
    types = np.unique(data['type'])
    for ty in types:
        data = Alldata[Alldata['type']==ty]
        t = data['dt']
        if ty == 'optical':
            y = data['mag']
            yerr = data['magerr']
            ylabel = 'mag'
        elif ty == 'x-ray' or ty == 'radio':
            y = data['flux']
            yerr = data['fluxerr']
            if ty == 'x-ray':
                ylabel = r'flux $(erg\cdot s^{-1}\cdot cm^{-2})$'
            else:
                ylabel = r'$\mu Jy$'
        else:
            raise KeyError('Invalid observation type!')
            
    

