import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import astropy.constants as const
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline

    

def HIloss_trapz(dwarf, Rlimit, check):

    testdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    testname = np.array(testdata['NAME'])

    if check == 'pearson':
        data = Table.read(dwarf+'.txt', format = 'ascii')
        radius = np.array(data['col1'])
        gas = np.array(np.log10(data['col2']))
    else:
        data = Table.read('data/'+dwarf+'_clean.dat', format = 'ascii')
        radius = np.array(data['R(kpc)'])
        gas = np.array(data['N_HI'])

    check = dwarf == testname
    test = testname[check]

    #This accounts for the 1.36 multiplicative factor based on He abundance
    if len(test) == 0:
        gas = gas + 0.13354
    else:
        gas = gas

    clean_radius = radius[gas < 50]
    clean_gas = gas[gas < 50]

    radiusfinal = np.linspace(min(clean_radius), max(clean_radius), len(clean_radius)*10)
    gasinterp = IUSpline(clean_radius, clean_gas, k=3)
    gasfinal = gasinterp(radiusfinal)

    if check == 'pearson':
        func = 0.5*np.pi*(radiusfinal*((3.08E21)**2))*(10**gasfinal)
    else:
        func = 2*np.pi*(radiusfinal*((3.08E21)**2))*(10**gasfinal)
        
    mass = np.trapz(func,radiusfinal)

    striplimit = Rlimit #np.max([Rlimit,np.min(radiusfinal[gasfinal < gaslimit])])
    stripradius = radiusfinal[radiusfinal < striplimit]
    stripfunc = func[radiusfinal < striplimit]
    nonstripmass = np.trapz(stripfunc, stripradius)

    frac = (mass - nonstripmass) / mass
    print dwarf+' = '+np.str(frac)

    gasmass = nonstripmass*(const.m_p / const.M_sun)

    return frac, gasmass, mass*(const.m_p / const.M_sun)

        

        

        

    
