import xspec as xs

def get_ctrt_to_flux(source_spec, energy_l, energy_h, nH, PhoIndex, get_unabs):
    xs.Xset.abund = 'wilm'
    xs.Fit.statMethod = 'cstat'
    spec = xs.Spectrum(source_spec)
    spec.ignore('**-%.1f %.1f-**'%(energy_l, energy_h))
    spec.background = None
    model = xs.Model('tbabs*cflux*powerlaw')
    if not get_unabs:
        model = xs.Model('cflux*tbabs*powerlaw')
    model.TBabs.nH.values = nH
    model.cflux.Emin = energy_l
    model.cflux.Emax = energy_h
    model.cflux.lg10Flux = -6.0
    model.powerlaw.PhoIndex = PhoIndex
    model.powerlaw.norm = 1.0
    model.powerlaw.norm.frozen = True
    fakeit_kwargs = {}
    fakeit_kwargs['exposure'] = 10000
    fakeit_kwargs['correction'] = 1.0
    fakeit_kwargs['backExposure'] = 1.0
    fakeit_kwargs['fileName'] = 'temp_fake.pha'
    xs.AllData.fakeit(1, xs.FakeitSettings(**fakeit_kwargs))
    spec = xs.AllData(1)
    spec.ignore('**-%.1f %.1f-**'%(energy_l, energy_h))
    spec.show()
    model.show()
    ctrt = spec.rate[-1]
    flux = 10**(model.cflux.lg10Flux.values[0])
    xs.AllModels.clear()
    xs.AllData.clear()
    return flux/ctrt


l, h = 4, 10
convert_ratio = get_ctrt_to_flux(source_spec='PC.pi',
                                 energy_l=l,
                                 energy_h=h,
                                 nH=1e-2,
                                 PhoIndex=1.7,
                                 get_unabs=True)
print('Energy ranging ({:.1f}-{:.1f}), factor = {}'.format(l,h,convert_ratio))
