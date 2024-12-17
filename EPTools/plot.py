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
    
    fig = plt.figure(figsize=(11,9),dpi=200)
    gs = fig.add_gridspec(len(types), hspace=0,height_ratios=np.array([1.5,2,1]))
    ax = gs.subplots(sharex=True)
    
    for i,ty in enumerate(types):
        data = Alldata[Alldata['type']==ty]
        t = data['dt']
        if ty == 'optical':
            y = data['mag']
            yerr = data['magerr']
            ylabel = 'mag'
            ax[i].invert_yaxis()
        elif ty == 'x-ray' or ty == 'radio':
            y = data['flux']
            yerr = data['fluxerr']
            if ty == 'x-ray':
                ylabel = r'flux $(erg\cdot s^{-1}\cdot cm^{-2})$'
            else:
                ylabel = r'$\mu Jy$'
        else:
            raise KeyError('Invalid observation type!')
        
        ax[i].errorbar(x=t,y=y,yerr=yerr)
        ax[i].set_ylabel(ylabel)
        if i == len(types)-1:
            ax[i].set_xlabel('t(day)')
            
    plt.savefig(output,dpi=300)
    

