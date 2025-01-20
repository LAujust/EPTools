from .utils import *
import dynesty
from fuzzywuzzy import fuzz

def xspec_fitting(sname,mname:str,grp=False,arf=None,rmf=None,rebin=5,instrument:str='WXT',plotmode='data',**fixed_par):
    """
    Args:
        sname[str/list]: spectrum file
        mname (str): model name
        grp (bool, optional): whether the file are grouped
        arf (str/list, optional): anxilury file
        rmf (str/list, optional): reponse file
        rebin (float, optional): rebin
        instrument (str, optional): 'WXT' or 'FXT'
        plotmode (str, optional): xspec plot mode (e.g. data, ldata, edata)

    Returns:
        tuple or list of tuple: (energies,edeltas,rates,errors,model,labels)
    """
    
    
    if instrument == 'WXT':
        el, eh = 0.5, 4
    elif instrument == 'FXT':
        el, eh = 0.5, 10
    else:
        raise KeyError('Input Valid Instrument Type')
    

    if isinstance(sname,str):
        if not grp:
            if not arf:
                arf = re.sub(r'\.(pha|pi)$', '.arf', sname)
            if not rmf:
                rmf = re.sub(r'\.(pha|pi)$', '.rmf', sname)
            s = xs.Spectrum(sname,arfFile=arf,respFile=rmf)
        else:
            s = xs.Spectrum(sname)
    elif isinstance(sname,list):
        for i, sn in enumerate(sname):
            if not grp:
                if not arf:
                    arf = re.sub(r'\.(pha|pi)$', '.arf', sname)
                if not rmf:
                    rmf = re.sub(r'\.(pha|pi)$', '.rmf', sname)
                s = xs.Spectrum(sname,arfFile=arf[i],respFile=rmf[i])
            else:
                s = xs.Spectrum(sname)
            
    xs.AllData.ignore("0.0-{:.1f} {:.1f}-**".format(el,eh))
    xs.AllData.ignore("bad")
            
    #Load Model
    m = xs.Model(mname)
    for key,value in fixed_par.items():
        exec('m.{}={}'.format(key,value))
        exec('m.{}.frozen=True'.format(key))
        
    xs.AllModels.show()
    xs.AllData.show()
    
    
    #Fit
    xs.Fit.renorm('auto')
    xs.Fit.nIterations = 300
    xs.Fit.statMethod = "cstat"
    xs.Fit.statMethod = "chi"
    xs.Fit.perform()
    xs.Fit.goodness(1000)
    xs.AllModels.calcFlux('{:.1f} {:.1f}'.format(el,eh))
    
    
    #Plot data
    xs.Plot.device = "/null"
    xs.Plot.xAxis="keV"
    xs.Plot.setRebin(rebin,20)
    xs.Plot(plotmode)
    
    if isinstance(sname,str):
        energies = xs.Plot.x(1,1)
        edeltas = xs.Plot.xErr(1,1)
        rates = xs.Plot.y(1,1)
        errors = xs.Plot.yErr(1,1)
        labels = xs.Plot.labels()
        model = xs.Plot.model()
        output = (energies,edeltas,rates,errors,model,labels)
    elif isinstance(sname,list):
        output = []
        for i in range(len(sname)):
            energies = xs.Plot.x(i+1,1)
            edeltas = xs.Plot.xErr(i+1,1)
            rates = xs.Plot.y(i+1,1)
            errors = xs.Plot.yErr(i+1,1)
            labels = xs.Plot.labels(i+1)
            model = xs.Plot.model(i+1)
            output.append((energies,edeltas,rates,errors,model,labels))
    
    xs.AllModels.clear()
    xs.AllData.clear()
    return output


def dynesty_fitting():
    pass

