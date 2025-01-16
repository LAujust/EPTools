from .utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl


def plot_gcn_data(file_dir='',output='standard_gcn.pdf',ttype=0):
    '''
    Plot GCN Circilar data and EP-WXT/FXT data from .csv file
    ttype[int]:    0 for days; 1 for seconds in logsapce

    Data Format:
    
    data.csv: [dt, terr, instrument, mag, magerr, flux, fluxerr, detection, type] 
    #flux in radio data often refers to flux density [mJy/uJy]
    #type should be either optical, x-ray or radio
    '''

    Alldata = pd.read_csv(file_dir).sort_values(by='dt').reset_index(drop=True)
    raw_types = np.unique(Alldata['type'])
    
    fig = plt.figure(figsize=(11,9),dpi=200)
    gs = fig.add_gridspec(len(raw_types), hspace=0)
    ax = gs.subplots(sharex=True)
    tscale = 24*3600 if ttype ==1 else 1

    ins_color = {'EP-WXT':'k','EP-FXT':'gray','XRT':'darkgreen',
                 'ATCA-5.5 GHz':'brown','ATCA-9 GHz':'indianred','e-MERLIN-5 GHz':'peru'}
    
    #Re-order types in an order of X-ray, Optical, Radio
    tyorder = {'x-ray': 0, 'optical': 1, 'radio': 2}
    types = sorted(raw_types, key=lambda element: tyorder[element])
    
    for i,ty in enumerate(types):
        data = Alldata[Alldata['type']==ty]
        t = data['dt']*tscale
        ax[i].grid()
        if ty == 'optical':
            det = np.array(data['detection'])
            uplim = [True if det[i]==0 else False for i in range(len(det))]
            bands = np.unique(data['band'])
            band = data['band']
            cmap = plt.cm.get_cmap('gist_rainbow', len(bands))
            band_color = {band:mpl.colors.rgb2hex(cmap(i)) for i,band in enumerate(bands)}
            color = [band_color[b] for b in band]
            y = data['mag']
            yerr = data['magerr']
            ylabel = 'mag'
            xerr = np.zeros(y.shape)
            ax[i].invert_yaxis()
            label = band

        elif ty == 'x-ray' or ty == 'radio':
            det = np.array(data['detection'])
            uplim = [True if det[i]==0 else False for i in range(len(det))]

            y = data['flux']
            yerr = data['fluxerr']
            xerr = data['terr']*tscale
            instrument = data['instrument']
            color = [ins_color[ins] for ins in instrument]
            label = instrument
            ax[i].set_yscale('log')
            if ty == 'x-ray':
                ylabel = r'flux $(erg\cdot s^{-1}\cdot cm^{-2})$'
            else:
                ylabel = r'$\mu Jy$'
        else:
            raise KeyError('Invalid observation type!')
        
        for ti,yi,yerri,xerri,colori,uplimi,labeli in zip(t,y,yerr,xerr,color,uplim,label):
            if uplimi:
                fmt = 'v'
                xerri, yerri = None, None
            else:
                fmt = 'o'
            ax[i].errorbar(x=ti,y=yi,yerr=yerri,xerr=xerri,fmt=fmt,markersize=5,capsize=2,color=colori,lolims=uplimi,label=labeli)
        ax[i].set_ylabel(ylabel)
        
        'Draw Legend'
        'https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend'
        handles, labels = ax[i].get_legend_handles_labels()
        handle_list, label_list = [], []
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        ncols = 2 if len(handle_list)>5 else 1
        ax[i].legend(handle_list, label_list, ncols=ncols)
        ax[i].autoscale()

        if i == len(types)-1:
            if ttype == 0:
                ax[i].set_xlabel('t(day)')
            elif ttype == 1:
                ax[i].set_xlabel('t(s)')
                ax[i].set_xscale('log')

    plt.suptitle(list(output.split('.'))[-2],fontsize=20)
    plt.savefig(output,dpi=300)
    

