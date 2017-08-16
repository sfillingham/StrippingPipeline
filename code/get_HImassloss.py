import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from scipy.integrate import cumtrapz
    

def HIloss(dwarf, Rlimit, gaslimit):

    data = Table.read(dwarf+'_clean.dat', format = 'ascii')
    radius = np.array(data['R(kpc)'])
    gas = np.array(data['N_HI'])

    clean_radius = radius[gas < 50]
    clean_gas = gas[gas < 50]

    lower = 0.0
    upper = clean_radius[0]
    mass = 0.0

    for i in range(len(clean_gas)-1):

        dr = (upper - lower)*(3.08E21) #Change units from kpc to cm
        r = clean_radius[i]*(3.08E21) #Change units from kpc to cm
        dm = 2*np.pi*r*(10**clean_gas[i])*dr

        mass = mass + dm

        lower = clean_radius[i]
        upper = clean_radius[i+1]
        #print clean_radius[i]


    striplimit = np.max([Rlimit,np.min(clean_radius[clean_gas < gaslimit])])
    print striplimit
    striplower = 0.0
    stripupper = clean_radius[0]
    stripmass = 0.0

    for i in range(len(clean_gas)-1):
        
        if clean_radius[i] < striplimit:

            dr = (stripupper - striplower)*(3.08E21) #Change units from kpc to cm
            r = clean_radius[i]*(3.08E21) #Change units from kpc to cm
            dm = 2*np.pi*r*(10**clean_gas[i])*dr

            stripmass = stripmass + dm

            lower = clean_radius[i]
            upper = clean_radius[i+1]

        else:
            stripmass = stripmass

    mass1 = cumtrapz(2*np.pi*clean_radius*((3.08E21)**2)*(10**clean_gas))
    print mass1[-1]

    cond = clean_radius == striplimit
    stripindex = np.where(cond)
    print stripindex
    
    print mass1[stripindex]

    #stripmass1 = cumtrapz(2*np.pi*r*(10**clean_gas[i])
    frac = stripmass / mass
    print mass
    print stripmass
    
    return frac

        

        

        

    
