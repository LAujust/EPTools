from .utils import *
import dynesty
import xspec as xs


def grp_data(sname,outputname,arf=None,rmf=None,group=1):
    cmd = ['grppha %s %s'%(sname,outputname)]
    comm = "comm='group min %s "%group
    if arf:
        comm += " & chkey ANCRFILE %s"%arf
    if rmf:
        comm += " & chkey RESPFILE %s"%rmf
    comm += " & bad 0-29 & exit'"
    cmd.append(comm)
    print(cmd)
    subprocess.run(cmd,capture_output=True, text=True)

def xspec_fitting(sname,mname:str,grp=False,arf=None,rmf=None,rebin=5,stat='cstat',instrument='WXT',untied=None,plotmode='data resid',chdir=None,**fixed_par):
    """
    !!!Single Spectrum Fitting or Simutaneously Fitting!!!
    !!!To fit single Spectrum, same should be a str; for simutaneously fitting, sname should be a 
    list of sname
    Args:
        sname[str/list]: spectrum file
        mname (str): model name
        grp (bool, optional): whether the file are grouped
        arf (str/list, optional): anxilury file
        rmf (str/list, optional): reponse file
        rebin (float, optional): rebin
        instrument (str or list, optional): 'WXT' or 'FXT'
        plotmode (str, optional): xspec plot mode (e.g. data, ldata, edata)
        untied (dict): untied parameters. The elements are parName:[mNum,value], e.g. {'tbabs.nH':[2,0.1]}
        chdir (str): change home dir
        fixed_par: fixed parameters. Second and third, forth model parameters are the same with the first one by default

    Returns:
        tuple or list of tuple: (energies,edeltas,rates,errors,model,labels)
    """
    if chdir:
        os.chdir(chdir)
    
    xs.Xset.allowNewAttributes = True

    erange = {'WXT':[0.5,4],'FXT':[0.5,10]}
    

    if isinstance(sname,str):
        if not grp:
            if not arf:
                arf = re.sub(r'\.(pha|pi)$', '.arf', sname)
            if not rmf:
                rmf = re.sub(r'\.(pha|pi)$', '.rmf', sname)
            s = xs.Spectrum(sname,arfFile=arf,respFile=rmf)
        else:
            s = xs.Spectrum(sname)
        xs.AllData.ignore("bad")
        xs.AllData.ignore("0.0-{:.1f} {:.1f}-**".format(erange[instrument][0],erange[instrument][1])) 
    elif isinstance(sname,list):
        for i, sn in enumerate(sname):
            if not grp:
                if not arf:
                    arf = re.sub(r'\.(pha|pi)$', '.arf', sname)
                if not rmf:
                    rmf = re.sub(r'\.(pha|pi)$', '.rmf', sname)

                si = xs.AllData(sn,arfFile=arf[i],respFile=rmf[i])
                si.ignore("bad")
                si.ignore("0.0-{:.1f} {:.1f}-**".format(erange[instrument[i]][0],erange[instrument[i]][1]))
                # xs.AllData(sn)
                # xs.AllData(i+1).response.arf = arf[i]
                # xs.AllData(i+1).response.rmf = rmf[i]
            else:
                xs.AllData('{}:{} {}'.format(i+1,i+1,sn))
                xs.AllData(i+1).ignore("bad")
                xs.AllData(i+1).ignore("0.0-{:.1f} {:.1f}-**".format(erange[instrument[i]][0],erange[instrument[i]][1]))
        
            
    #Load Model and freeze parameters
    m = xs.Model(mname)
    # if isinstance(sname,list):
    #     for i in range(len(sname)):
    #         exec('m%s = AllModels(%s)'%(i+1,i+1))
    for key,value in fixed_par.items():
        exec('m.{}={}'.format(key,value))
        exec('m.{}.frozen=True'.format(key))
    if untied:
        for par,value in untied.items():
            exec('xs.AllModels({}).{}={}'.format(value[0],par,value[1]))
        
    xs.AllModels.show()
    xs.AllData.show()
    
    #Fit
    xs.Fit.renorm('auto')
    xs.Fit.nIterations = 300
    xs.Fit.statMethod = stat
    xs.Fit.perform()
    # xs.Fit.goodness(200)
    if isinstance(sname,str):
        xs.AllModels.calcFlux('{:.1f} {:.1f}'.format(erange[instrument][0],erange[instrument][1]))
    else:
        for i in range(len(sname)):
            xs.AllModels.calcFlux('{:.1f} {:.1f}'.format(erange[instrument[i]][0],erange[instrument[i]][1]))
    
    
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
        resid = xs.Plot.y(1,2)
        residerr = xs.Plot.yErr(1,2)
        output = (energies,edeltas,rates,errors,model,resid,residerr,labels)
    elif isinstance(sname,list):
        output = []
        labels = xs.Plot.labels()
        for i in range(len(sname)):
            energies = xs.Plot.x(i+1,1)
            edeltas = xs.Plot.xErr(i+1,1)
            rates = xs.Plot.y(i+1,1)
            errors = xs.Plot.yErr(i+1,1)
            model = xs.Plot.model(i+1)
            resid = xs.Plot.y(i+1,2)
            residerr = xs.Plot.yErr(i+1,2)
            output.append((energies,edeltas,rates,errors,model,resid,residerr,labels))
    
    xs.AllModels.clear()
    xs.AllData.clear()
    return output


def dynesty_fitting(prior_transform, log_likeli,bounds,nlive,log_args):

    pass

