from .utils import *

__all__ = ['grp_data','xspec_fitting','dynesty_fitting']

MODEL_COMP_PARAM = {
    'TBabs':['nH'],
    'cflux':['Emin','Emax','lg10Flux'],
    'powerlaw':['PhoIndex','norm'],
    'bbody':['kT','norm'],
    'apec':['Abundanc','kT','Redshift','norm']
}


def grp_data(src,bkg,outputname='PC.pi',arf=None,rmf=None,group=1):
    os.system(f"grppha infile={src} outfile=PC.pi chatter=0 comm='group min {group} & chkey RESPFILE {rmf} & chkey ANCRFILE {arf} & chkey BACKFILE {bkg} & exit'")

def xspec_fitting(sname,mname:str,grp=False,arf=None,rmf=None,rebin=None,stat='cstat',instrument='WXT',untied=None,plotmode='data resid',chdir=None,N=500,chatter=10,**fixed_par):
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
        untied (dict): untied parameters. e.g. ['powerlaw.norm=1,0.1',...]
        chdir (str): change home dir
        fixed_par: fixed parameters. Second and third, forth model parameters are the same with the first one by default

    Returns:
        tuple or list of tuple: (energies,edeltas,rates,errors,model,labels)
    """
    if chdir:
        os.chdir(chdir)
    
    xs.Xset.allowNewAttributes = True
    xs.Xset.chatter = chatter

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
        
    xs.AllData.ignore("bad")
            
    #Load Model and freeze parameters
    m = xs.Model(mname)
    # if isinstance(sname,list):
    #     for i in range(len(sname)):
    #         exec('m%s = AllModels(%s)'%(i+1,i+1))
    for key,value in fixed_par.items():
        exec('m.{}={}'.format(key,value))
        exec('m.{}.frozen=True'.format(key))
    if untied:
        for item in untied:
            key, value = item.split("=")
            values = [float(v) if "." in v else int(v) for v in value.split(",")]  # 转换成数字
            exec('xs.AllModels({}).{}={}'.format(values[0],key,values[1]))
        
    xs.AllModels.show()
    xs.AllData.show()
    #xs.Xset.chatter = 10
    
    #Fit
    xs.Fit.renorm('auto')
    xs.Fit.nIterations = N
    xs.Fit.statMethod = stat
    xs.Fit.perform()
    xs.Fit.show()
    cstat = xs.Fit.statistic
    dof = xs.Fit.dof
    reduced_cstat = cstat / dof
    # xs.Fit.goodness(200)

    # Split on *, (, ), +
    parts = re.split(r'[\*\(\)\+]', mname)
    comps = [p for p in parts if p]

    par_table = Table()
    par_table['name'] = [mname]
    par_table['sname'] = [sname]
    par_table['cstat/dof'] = reduced_cstat

    # comps = m.componentNames
    # print(comps)
    # for comp in comps:
    #     pars = []
    #     print(comp)
    #     exec("pars = m.%s.parameterNames"%comp)
    #     print(pars)
    #     for par in pars:
    #         exec("par_table[par]=[m.%s.%s.values[0]]"%(comp,par))

    par_num = 0
    for comp in comps:
        if comp == 'tbabs':
            comp = 'TBabs'
        pars = MODEL_COMP_PARAM[comp]
        for par in pars:
            par_num += 1
            xs.Fit.error(str(par_num))
            #print("print(par,m.%s.%s.error)"%(comp,par))
            #exec("print(par,m.%s.%s.values)"%(comp,par))
            #exec("print(par,m.%s.%s.error)"%(comp,par))
            exec("par_table[par]=[m.%s.%s.values[0]]"%(comp,par))
            exec("par_table['%s_err_low']=[m.%s.%s.error[0]]"%(par,comp,par))
            exec("par_table['%s_err_high']=[m.%s.%s.error[1]]"%(par,comp,par))
    
    par_table.pprint(max_lines=-1,max_width=-1)

    if isinstance(sname,str):
        xs.AllModels.calcFlux('{:.1f} {:.1f}'.format(erange[instrument][0],erange[instrument][1]))
        xs.AllModels.calcFlux('0.5 2')
        xs.AllModels.calcFlux('0.5 {:.1f}'.format(erange[instrument][1]))

    else:
        for i in range(len(sname)):
            xs.AllModels.calcFlux('{:.1f} {:.1f}'.format(erange[instrument[i]][0],erange[instrument[i]][1]))
            xs.AllModels.calcFlux('0.5 2')
            xs.AllModels.calcFlux('0.5 {:.1f}'.format(erange[instrument][i][1]))

    
    
    #Plot data
    xs.Plot.device = "/null"
    xs.Plot.xAxis="keV"
    if rebin:
        if len(rebin)>1:
            xs.Plot.setRebin(minSig=rebin[0],maxBins=rebin[1])
        else:
            xs.Plot.setRebin(rebin[0])
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
    return output, par_table


def dynesty_fitting(prior_transform, log_likeli,bounds,nlive,log_args):

    pass

