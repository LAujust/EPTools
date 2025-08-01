from .utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits

DEFAULT_COLORS = ['#002FA7','#FFE76F',
                  '#01847F','#F9D2E4',
                  '#FF0000','#492D22',
                  '#FF770F','#000026',
                  '#003153','#7356B1',
                  '#BDB5D7','#4C7543',
                  '#194955','#F49227']


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
    



def xspec_plot(data,save_dir=None,leg=None,color='random',plotstyle='step'):

    if isinstance(data,tuple):
        energies,edeltas,rates,errors,model,resid,residerr,labels = data
        nE = len(energies)
        stepenergies = list()
        for i in range(nE):
            stepenergies.append(energies[i] - edeltas[i])
        stepenergies.append(energies[-1]+edeltas[-1])
        model.append(model[-1])
        #resid.append(resid[-1])

        #plt.yscale('log')
        fig = plt.figure(figsize=(7,6))
        gs = fig.add_gridspec(2, hspace=0, height_ratios=[3,1])
        ax = gs.subplots(sharex=True)
        ax[0].set_ylabel(labels[1])
        ax[1].set_xlabel(labels[0])
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[1].set_ylabel('Residual')
        ax[0].set_title(labels[2])
        ax[0].errorbar(energies,rates,xerr=edeltas,yerr=np.abs(errors),fmt='.',color='dimgrey',label=leg)
        ax[1].errorbar(energies,resid,xerr=edeltas,yerr=np.abs(residerr),color='dimgrey',fmt='.')
        if plotstyle == 'step':
            ax[0].step(stepenergies,model,where='post',color='royalblue')
        elif plotstyle == 'line':
            ax[0].plot(stepenergies,model,color='royalblue')
        else:
            raise KeyError('Not valid plotstyle (step/line)')
        
        ax[0].legend()
        plt.grid()
        plt.tight_layout()
        if save_dir:
            plt.savefig(save_dir,dpi=300)
        else:
            plt.show()


    elif isinstance(data,list) and isinstance(data[0],tuple):
        rows = 2
        fig = plt.figure(figsize=(7,6))
        gs = fig.add_gridspec(rows, hspace=0, height_ratios=[3,1])
        ax = gs.subplots(sharex=True)
        if color == 'random':
            cmap = plt.cm.get_cmap('gist_rainbow',30)
            rand_color_idx = np.random.randint(0,30,len(data))
            random_color = [mpl.colors.rgb2hex(cmap(i)) for i in rand_color_idx]
        elif color == 'default':
            random_color = np.random.choice(DEFAULT_COLORS,len(data))

        for i in range(len(data)):
            energies,edeltas,rates,errors,model,resid,residerr,labels = data[i]
            nE = len(energies)
            stepenergies = list()
            for j in range(nE):
                stepenergies.append(energies[j] - edeltas[j])
            stepenergies.append(energies[-1]+edeltas[-1])
            model.append(model[-1])

            ax[0].errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.',color=random_color[i],label=leg[i])
            if plotstyle == 'step':
                ax[0].step(stepenergies,model,where='post',color='dimgrey')
            elif plotstyle == 'line':
                ax[0].plot(stepenergies,model,color='royalblue')
            else:
                raise KeyError('Not valid plotstyle (step/line)')
            ax[1].errorbar(energies,resid,xerr=edeltas,yerr=residerr,color=random_color[i],fmt='.')
            ax[0].legend()
            ax[0].set_ylabel(labels[1])
            ax[1].set_xlabel(labels[0])
            ax[0].set_xscale('log')
            ax[0].set_yscale('log')
            ax[1].set_ylabel('Residual')
            ax[0].set_title(labels[2])
            
        #ax[1].hlines(0.0,0.0,10.0,ls='dashed',color='maroon')
        ax[0].grid()
        ax[1].grid()

        if save_dir:
            plt.savefig(save_dir,dpi=300)
        else:
            plt.show()


def lcurve_plot(src,bkg,save_dir=None,binsize=10,scale=1./12,rx=None,sep=False,show=False):
    with fits.open(src) as hdu:
        TSTART = hdu[0].header['TSTART']
        T0 = TSTART
        #DATE_OBS = hdu[0].header['DATE-OBS']
        data = hdu[1].data
        TIME = data['TIME']
        RATE = data['RATE']
        ERROR = data['ERROR']
    with fits.open(bkg) as hdu:
        bkg = hdu[1].data
        TIME_bkg = bkg['TIME']
        RATE_bkg = bkg['RATE']
        ERROR_bkg = bkg['ERROR']

    t,rate,error = [],[],[]
    t_bkg,rate_bkg,error_bkg = [],[],[] #Scaled
    #Rebin
    pin = TIME[0]
    while pin+binsize < TIME[-1]:
        idx = np.where((TIME>pin) & (TIME<pin+binsize))[0]
        true_size = len(idx)
        if true_size == 0:
            pin += binsize
            continue
        else:
            t.append(sum(TIME[idx])/true_size)
            rate.append(sum(RATE[idx])/true_size)
            error.append(np.sqrt(sum(ERROR[idx]**2))/true_size)
            rate_bkg.append(scale*sum(RATE_bkg[idx])/true_size)
            error_bkg.append(scale*np.sqrt(sum(ERROR_bkg[idx]**2))/true_size)
            pin += binsize
            
            
    if sep:
        fig = plt.figure()
        gs = fig.add_gridspec(3,hspace=0)
        ax = gs.subplots(sharex=True)
        ax[0].errorbar(t,rate,yerr=error,xerr=binsize/2,color='steelblue',fmt='.',alpha=0.7,label='Src')
        ax[2].errorbar(t,rate_bkg,yerr=error_bkg,xerr=binsize/2,color='grey',alpha=0.7,fmt='.',label='Scaled bkg')
        ax[1].errorbar(t,np.array(rate)-np.array(rate_bkg),yerr=np.sqrt(np.array(error)**2+np.array(error_bkg)**2),
                    xerr=binsize/2,
                    color='darkorange',alpha=0.7,fmt='.',label='Net')
        if rx:
                ax[2].set_xlim(rx)
        
        for i in range(3):
            ax[i].hlines(0,t[0],t[-1],color='k',ls='--')
            ax[i].legend()
            ax[i].grid()
            ax[i].set_ylabel('counts/s')

        ax[2].set_xlabel('$\mathrm{T-T_{0}}=$'+'{} (bintime={:.1f}s)'.format(T0,binsize))
        if save_dir:
            plt.savefig(save_dir,dpi=300)
        if show:
            plt.show()
        else:
            return fig, ax
    
    else:
        fig, ax = plt.subplots(dpi=100,figsize=(7,5))
        # gs = fig.add_gridspec(2, hspace=0,height_ratios=[1.5,1])
        # ax = gs.subplots(sharex=True)

        ax.errorbar(t,rate,yerr=error,xerr=binsize/2,color='steelblue',fmt='.',alpha=0.7,label='Src')
        ax.errorbar(t,rate_bkg,yerr=error_bkg,xerr=binsize/2,color='grey',alpha=0.7,fmt='.',label='Scaled bkg')
        ax.errorbar(t,np.array(rate)-np.array(rate_bkg),yerr=np.sqrt(np.array(error)**2+np.array(error_bkg)**2),
                    xerr=binsize/2,
                    color='darkorange',alpha=0.7,fmt='.',label='Net')
        ax.hlines(0,t[0],t[-1],color='k',ls='--')
        ax.legend()
        ax.grid()
        ax.set_ylabel('counts/s')
        ax.set_xlabel('$\mathrm{T-T_{0}}=$'+'{} (bintime={:.1f}s)'.format(T0,binsize))
        if rx:
            ax.set_xlim(rx)
        if save_dir:
            plt.savefig(save_dir,dpi=300)
        if show:
            plt.show()
        else:
            return fig, ax

