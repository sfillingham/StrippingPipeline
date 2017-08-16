
import numpy as np
from astropy.table import Table
import astropy.units as u
import astropy.constants as consts

def LT(dwarf):
    
    xx = np.str(dwarf)
    datafile = open(xx+'.dat','r')
    data = datafile.readlines()
    
    N = np.size(data)
    radius = np.ndarray(N-2)
    columndensity = np.ndarray(N-2)

    for t in range(N-2):
        radius[t] = float(data[t+2].strip().split()[1])
        columndensity[t] = float(data[t+2].strip().split()[3])


    output = Table([radius, columndensity], names=('R(kpc)', 'N_HI'), meta={'HI Surface Density': 'LittleThings'})

    outputfile = xx+'_clean.dat'
    output.write(outputfile, format = 'ascii')
    
    return radius, columndensity

def THINGS():

    data = Table.read('data/THINGS_clean.dat', format = 'ascii')
    altdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    dname = np.array(data['NAME'])
    r = np.array(data['R(kpc)'])
    sigma = np.array(data['SigmaHI'])
    altname = np.array(altdata['altname'])
    outname = np.array(altdata['NAME'])
    
    print dname
    print altname
    

    for i in range(len(altname)):

        dwarf = altname[i]
        cut = dname == dwarf
        radius = r[cut]
        sig = np.log10(sigma[cut]*(consts.M_sun/consts.m_p)*(1e-2/consts.pc)**2*u.m**2) #final units are log(atoms cm^-2)
        

        output = Table([radius, sig], names=('R(kpc)', 'N_HI'), meta={'HI Surface Density': 'LittleThings'})
    
        outputfile = 'data/'+outname[i]+'_clean.dat'
        output.write(outputfile, format = 'ascii')


def SHIELD():

    data = Table.read('SHIELD_list.dat', format = 'ascii')
    dname = np.array(data['NAME'])
    
    for i in range(len(dname)):

        objdata = Table.read('data/'+dname[i]+'.dat', format = 'ascii')
        radius = np.array(objdata['R(kpc)'])
        sigma = 10**np.array(objdata['N_HI'])

        sig = np.log10(sigma*(consts.M_sun/consts.m_p)*(1e-2/consts.kpc)**2*u.m**2) #final units are log(atoms cm^-2)
        

        output = Table([radius, sig], names=('R(kpc)', 'N_HI'), meta={'HI Surface Density': 'SHIELD'})
    
        outputfile = 'data/'+dname[i]+'_clean_kpc.dat'
        output.write(outputfile, format = 'ascii')







