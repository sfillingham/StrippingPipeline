import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import astropy.units as u
import astropy.constants as consts
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline
import massprofile as mass

def strip(dwarf, massinput, concentration, r0, rho0, rs, rhos, hotdensity, velinfall, check):


    testdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    testname = np.array(testdata['NAME'])

    if check == 'pearson':
        data = Table.read(dwarf+'.txt', format = 'ascii')
        gasdata = np.array(np.log10(data['col2']))
        gas = gasdata[gasdata < 50]
        rdata = np.array(data['col1'])
        r = rdata[gasdata < 50]*(3.086e21) #change units from kpc to cm
        rfinal = np.linspace(min(r), max(r), len(r)*10)
    else:
        data = Table.read('data/'+dwarf+'_clean.dat', format = 'ascii')
        gasdata = np.array(data['N_HI'])
        gas = gasdata[gasdata < 50]
        rdata = np.array(data['R(kpc)'])
        r = rdata[gasdata < 50]*(3.086e21) #change units from kpc to cm
        rfinal = np.linspace(min(r), max(r), len(r)*10)

    jackass = dwarf == testname
    test = testname[jackass]

    #This accounts for the 1.36 multiplicative factor based on He abundance 
    if len(test) == 0:
        gas = gas + 0.13354
    else:
        gas = gas
    
    
    #print 'rfinal = '+np.str(rfinal/(3.086e21))
    gasinterp = IUSpline(r, gas, k=3)
    gasfinal = gasinterp(rfinal)
    
    strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
    
####################################################################################
## Determine the mass profile and restoring pressure, returning the stripping radius
####################################################################################
    
    if check == 'Apace':
    
    #This restoring pressure is based on a burkert profile for the dark matter halo by specifying both the virial mass and concentration.

        burkertmass = mass.burkert(r0, rho0, rfinal/(3.086e21))
        burkertmassrestore = (((10**gasfinal)*consts.G*burkertmass)/(rfinal**2))*(1e-4*consts.M_sun)*(u.s**2/u.m**3)
    
    #This restoring pressure is based on a burkert profile for the dark matter halo by specifying both the virial mass and concentration.

        NFWmass_pace = mass.NFW_pace(rs, rhos, rfinal/(3.086e21))
        NFWmassrestore_pace = (((10**gasfinal)*consts.G*NFWmass_pace)/(rfinal**2))*(1e-4*consts.M_sun)*(u.s**2/u.m**3)
        
        NFWcond = strip_pressure < NFWmassrestore_pace
        Bcond = strip_pressure < burkertmassrestore
    
        NFWradius = rfinal[NFWcond]/(3.086e21)
        Bradius = rfinal[Bcond]/(3.086e21)

        test1 = (len(NFWradius) > 0) & (len(Bradius) > 0)
        test2 = (len(NFWradius) == 0) & (len(Bradius) > 0)
        test3 = (len(NFWradius) > 0) & (len(Bradius) == 0)

        if test1:

            NFWcut = ((rfinal/(3.086e21)) == NFWradius[-1])
            NFWmass = NFWmass_pace[NFWcut]
            Bcut = ((rfinal/(3.086e21)) == Bradius[-1])
            Bmass = burkertmass[Bcut]
            return NFWradius[-1], Bradius[-1], NFWmass, Bmass
    
        elif test2:

            NFWcut = ((rfinal/(3.086e21)) == np.min(rfinal))
            NFWmass = NFWmass_pace[NFWcut]
            Bcut = ((rfinal/(3.086e21)) == Bradius[-1])
            Bmass = burkertmass[Bcut]
            return np.min(rfinal)/(3.086e21), Bradius[-1], NFWmass, Bmass

        elif test3:

            NFWcut = ((rfinal/(3.086e21)) == NFWradius[-1])
            NFWmass = NFWmass_pace[NFWcut]
            Bcut = ((rfinal/(3.086e21)) == np.min(rfinal))
            Bmass = burkertmass[Bcut]
            return NFWradius[-1], np.min(rfinal)/(3.086e21), NFWmass, Bmass

        else:

            NFWcut = ((rfinal/(3.086e21)) == np.min(rfinal))
            NFWmass = NFWmass_pace[NFWcut]
            Bcut = ((rfinal/(3.086e21)) == np.min(rfinal))
            Bmass = burkertmass[Bcut]
            return np.min(rfinal)/(3.086e21), np.min(rfinal)/(3.086e21), NFWmass, Bmass
    
    else:
    #This restoring pressure is based on a NFW profile for the dark matter halo by specifying both the virial mass and concentration.

        NFWmass = mass.NFW(massinput, concentration, rfinal/(3.086e21))
        NFWrestore = (((10**gasfinal)*consts.G*NFWmass)/(rfinal**2))*(1e-4*consts.M_sun)*(u.s**2/u.m**3)
        NFWcond = strip_pressure < NFWrestore
    
        DMradius = rfinal[NFWcond]
        NFWradius = rfinal[NFWcond]/(3.086e21)

        test1 = (len(NFWradius) > 0)

        if test1:

            DMcut = rfinal == DMradius[-1]
            DMmass = NFWmass[DMcut]
            return NFWradius[-1], DMmass#, 10**gasfinal[NFWcond], NFWradius, NFWmass[NFWcond]

        else:

            DMcut = rfinal == np.min(rfinal)
            DMmass = NFWmass[DMcut]
            return np.min(rfinal)/(3.086e21), DMmass

