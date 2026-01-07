from .utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.nddata import Cutout2D
from matplotlib.patches import Circle
import plotly.express as px
import re

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
    



def xspec_plot(data,save_dir=None,title=None,color='random',plotstyle='step',model_color='#4575b4',model_leg=None):

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
        if not title:
            ax[0].set_title(labels[2])
        else:
            ax[0].set_title(labels[2]+' (%s)'%title)
        ax[0].errorbar(energies,rates,xerr=edeltas,yerr=np.abs(errors),fmt='.',color='dimgrey')
        ax[1].errorbar(energies,resid,xerr=edeltas,yerr=np.abs(residerr),color='dimgrey',fmt='.')
        if plotstyle == 'step':
            ax[0].step(stepenergies,model,where='post',color=model_color,label=model_leg)
        elif plotstyle == 'line':
            ax[0].plot(stepenergies,model,color=model_color,label=model_leg)
        else:
            raise KeyError('Not valid plotstyle (step/line)')
        
        ax[0].legend()
        xmin, xmax = ax[0].get_xlim()
        ax[0].set_xlim([np.max([xmin,0.5]),np.min([xmax,10])])
        ax[0].set_ylim([1e-8,3e0])
        ax[0].tick_params(axis='both',which='both',direction='in')
        ax[1].tick_params(axis='both',which='both',direction='in')
        plt.grid()
        plt.tight_layout()
        if save_dir:
            plt.savefig(save_dir,bbox_inches='tight',dpi=300)
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

            ax[0].errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.',color=random_color[i])
            if plotstyle == 'step':
                ax[0].step(stepenergies,model,where='post',color=model_color,label=model_leg)
            elif plotstyle == 'line':
                ax[0].plot(stepenergies,model,color=model_color,label=model_leg)
            else:
                raise KeyError('Not valid plotstyle (step/line)')
            ax[1].errorbar(energies,resid,xerr=edeltas,yerr=residerr,color=random_color[i],fmt='.')
            ax[0].legend()
            ax[0].set_ylabel(labels[1])
            ax[1].set_xlabel(labels[0])
            ax[0].set_xscale('log')
            ax[0].set_yscale('log')
            ax[1].set_ylabel('Residual')
            ax[0].tick_params(axis='both',which='both',direction='in')
            ax[1].tick_params(axis='both',which='both',direction='in')
            if not title:
                ax[0].set_title(labels[2])
            else:
                ax[0].set_title(labels[2]+' (%s)'%title)
            
        #ax[1].hlines(0.0,0.0,10.0,ls='dashed',color='maroon')
        ax[0].grid()
        ax[1].grid()

        if save_dir:
            plt.savefig(save_dir,bbox_inches='tight',dpi=300)
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

        T0 = Time('2020-01-01 00:00:00') + TSTART*u.second
        T0 = T0.iso
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
            
            
    #correct t
    t = np.array(t)
    if np.min(t)>1e8:
        t -= TSTART
            
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
            ax[i].tick_params(axis='both',which='both',direction='in')

        ax[2].set_xlabel('$\mathrm{T-T_{0}}=$'+'{} (bintime={:.1f}s)'.format(T0,binsize))
        if save_dir:
            plt.savefig(save_dir,bbox_inches='tight',dpi=300)
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
        ax.tick_params(axis='both',which='both',direction='in')
        ax.set_ylabel('counts/s')
        ax.set_xlabel('$\mathrm{T-T_{0}}=$'+'{} (bintime={:.1f}s)'.format(T0,binsize))
        if rx:
            ax.set_xlim(rx)
        if save_dir:
            plt.savefig(save_dir,bbox_inches='tight',dpi=300)
        if show:
            plt.show()
        else:
            return fig, ax



from matplotlib.patches import Circle

def plot_fits_with_region(fits_dir, reg_file, contrast=0.2, plot_method='matplotlib',size=1, output_file=None):
    """Plot a FITS image with a specified region from a DS9 region file.

    Args:
        fits_dir (str): image FITS file path.
        reg_file (str): region file path in DS9 format.
        contrast (float, optional): zscale constrast. Defaults to 0.2.
        plot_method (str, optional): Defaults to 'matplotlib'. Options: 'matplotlib', 'plotly'.
        size (int, optional): image size to plot in degree. Defaults to 1.
        output_file (str, optional): output image. Defaults to None.

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    # Read the FITS file
    with fits.open(fits_dir) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

     # Read the region file
    with open(reg_file, 'r') as f:
        lines = [line for line in f.readlines() if line.strip().lower().startswith('circle')]
        print(lines)

        match = re.search(r'circle\s*\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^\)"]+)', lines[0])
        if match:
            ra, dec, radius = map(float, match.groups())
        else:
            raise ValueError("Region file format not recognized. Expected DS9 circle format.")

    # Create a SkyCoord object for the region center
    region_center = SkyCoord(ra, dec, unit='deg')

    # Perform a cutout of the image
    cutout_size = (size * u.deg, size * u.deg) 
    cutout = Cutout2D(data, region_center, cutout_size, wcs=wcs)

    # Normalize the image using ZScale
    zscale = ZScaleInterval(contrast=contrast)
    vmin, vmax = zscale.get_limits(cutout.data)

    # Plot using the specified method
    if plot_method == 'matplotlib':
        fig, ax = plt.subplots(figsize=(7,5),subplot_kw={'projection': cutout.wcs})
        im = ax.imshow(cutout.data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        ax.set_title('FITS Image with Region')
        # Add the region marker
        region_pix = cutout.wcs.world_to_pixel(region_center)
        radius_pix = radius / cutout.wcs.pixel_scale_matrix[0, 0]
        circle = Circle((region_center.ra.deg,region_center.dec.deg),radius/3600, edgecolor='red', facecolor='none', lw=1, transform=ax.get_transform('world'))
        ax.add_patch(circle)
        #plt.colorbar(im, ax=ax, label='ctr')
        if output_file:
            plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.show()

    elif plot_method == 'plotly':
        # Create a plotly figure
        fig = px.imshow(cutout.data, color_continuous_scale='gray', origin='lower', zmin=vmin, zmax=vmax)
        fig.update_layout(title='FITS Image with Region', xaxis_title='RA', yaxis_title='Dec')
        # Add the region marker
        region_pix = cutout.wcs.world_to_pixel(region_center)
        fig.add_shape(type='circle',
                      xref='x', yref='y',
                      x0=region_pix[0] - radius_pix, x1=region_pix[0] + radius_pix,
                      y0=region_pix[1] - radius_pix, y1=region_pix[1] + radius_pix,
                      line=dict(color='red'))
        if output_file:
            fig.write_html(output_file)
        fig.show()

    else:
        raise ValueError("Invalid plot method. Choose 'matplotlib' or 'plotly'.")