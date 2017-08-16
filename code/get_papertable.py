import numpy as np
from astropy.table import Table


def table():

    LTdata = Table.read('LT_galaxylist.dat', format = 'ascii')
    Tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    Pdata = Table.read('best_fit_parameters_mlfixed.dat', format = 'ascii')

    LTname = np.array(LTdata['NAME'])
    Tname = np.array(Tdata['NAME'])
    Mv = np.array(LTdata['Mv'])
    LTmstar = ((4.83-(Mv))/2.5)
    Tmstar = np.array(Tdata['Mstar'])
    r0 = np.array(Pdata['r0'])
    rs = np.array(Pdata['rs'])
    rho0 = np.array(Pdata['rho0'])
    rhos = np.array(Pdata['rhos'])
    LTc = np.array(LTdata['c'])
    LTmvir = np.array(LTdata['virmass'])

    
