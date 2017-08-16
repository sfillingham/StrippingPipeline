import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline

def strip(lowmass, highmass, lowvmax, highvmax):

    lowdata = Table.read(lowmass+'_clean.dat', format = 'ascii')
    highdata = Table.read(highmass+'_clean.dat', format = 'ascii')

    lowgasdata = np.array(lowdata['N_HI'])
    lowgas = lowgasdata[lowgasdata < 50]
    lowrdata = np.array(lowdata['R(kpc)'])
    lowr = lowrdata[lowgasdata < 50]*(3.086e21) #change units from kpc to cm
    #print len(lowr)
    lowrfinal = np.linspace(min(lowr), max(lowr), len(lowr)*10)
    lowgasinterp = IUSpline(lowr, lowgas, k=3)
    lowgasfinal = lowgasinterp(lowrfinal)

    highgasdata = np.array(highdata['N_HI'])
    highgas = highgasdata[highgasdata < 50]
    highrdata = np.array(highdata['R(kpc)'])
    highr = highrdata[highgasdata < 50]*(3.086e21) #change units from kpc to cm
    #print len(highr)
    highrfinal = np.linspace(min(highr), max(highr), len(highr)*10)
    highgasinterp = IUSpline(highr, highgas, k=3)
    highgasfinal = highgasinterp(highrfinal)

    strip_pressure = (10**(-3.5))*(4e4) #units: cm**(-3) (km/s)**2
    print strip_pressure
    protonmass = 1.673e-27 #kg

    lowrestore = (2.0*(lowvmax**2)*(10**lowgasfinal)) / (lowrfinal)
    #print lowrestore
    highrestore = (2.0*(highvmax**2)*(10**highgasfinal)) / (highrfinal)
    #print highrestore

    lowcond = strip_pressure > lowrestore
    highcond = strip_pressure > highrestore

    lowradius = lowrfinal[lowcond]
    highradius = highrfinal[highcond]

    plt.figure(1)
    plt.plot(highr, highgas, 'bo')
    plt.plot(highrfinal, highgasfinal, 'g-')

    plt.figure(2)
    plt.plot(lowr, lowgas, 'bo')
    plt.plot(lowrfinal, lowgasfinal, 'g-')

    return np.min(lowradius)/(3.086e21), np.min(highradius)/(3.086e21)

